from ensemble_analyser.conformer import Conformer
from ensemble_analyser.ioFile import read_ensemble, save_snapshot
from ensemble_analyser.logger import create_log, ordinal
from ensemble_analyser.parser_arguments import parser_arguments
from ensemble_analyser.parser_parameter import get_conf_parameters
from ensemble_analyser.IOsystem import SerialiseEncoder
from ensemble_analyser.protocol import Protocol, load_protocol
from ensemble_analyser.pruning import calculate_rel_energies, check_ensemble
from ensemble_analyser.grapher import Graph

import ase
import time
import json
import datetime
from tabulate import tabulate
import os


MAX_TRY = 5


def launch(
    idx, conf, protocol, cpu, log, temp, ensemble, try_num: int = 1
) -> None:
    """
    Run the calculation for each conformer

    idx | int : index of the calculation
    conf | Conformer : conformer instance
    protocol | Protocol : protocol instance
    cpu | int : number of cpu to allocate
    log : logger instance
    temp | float : temperature [K]
    ensemble | list : whole ensemble list

    return None
    """

    log.info(
        f"{idx}. Running {ordinal(int(protocol.number))} PROTOCOL -> CONF{conf.number}"
    )
    try:
        st = time.perf_counter()

        calculator, label = protocol.get_calculator(
            cpu=cpu, charge=conf.charge, mult=conf.mult
        )
        atm = conf.get_ase_atoms(calculator)
        try:
            atm.get_potential_energy()
        except ase.calculators.calculator.PropertyNotImplementedError:
            pass

        end = time.perf_counter()

        os.rename(
            f"{label}.out", f"{conf.folder}/protocol_{protocol.number}.out"
        )
        if protocol.freq:
            os.rename(
                f"{label}.hess", f"{conf.folder}/protocol_{protocol.number}.hess"
            )
        os.remove(f"{label}.gbw")

    except ase.calculators.calculator.CalculationFailed:
        with open(f"{label}.out") as f:
            fl = f.read()
        log.error("\n".join(fl.splitlines()[-6:-3]))
        log.critical(
            f"\n{'='*20}\nCRITICAL ERROR\n{'='*20}\nSome sort of error have been encountered during the calculation of the calculator.\n{'='*20}\nExiting\n{'='*20}\n"
        )
        raise RuntimeError(
            "Some sort of error have been encountered during the calculation of the calculator."
        )

    if not get_conf_parameters(
        conf, protocol.number, protocol, end - st, temp, log
    ):
        if try_num <= MAX_TRY:
            log.error(
                f"ERROR: During calculation of CONF_{conf.number} a server error occur and the energy could not be parsed; re-running protocol {protocol.number} on the same conformer for the {ordinal(try_num)} time"
            )
            time.sleep(10)
            launch(
                idx,
                conf,
                protocol,
                cpu,
                log,
                temp,
                ensemble,
                try_num=try_num + 1,
            )
        else:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nMax number of re-run ({MAX_TRY}) executed for CONF_{conf.number}.{'='*20}\nExiting\n{'='*20}"
            )
            raise RuntimeError(
                f"Max number of re-run ({MAX_TRY}) executed for CONF_{conf.number}. Exiting"
            )

    json.dump(
        {i.number: i.__dict__ for i in ensemble},
        open("checkpoint.json", "w"),
        indent=4,
        cls=SerialiseEncoder,
    )


def run_protocol(conformers, p, temperature, cpu, log) -> None:
    """
    Run the protocol for each conformer

    conformers | list : whole ensemble list
    p | Protocol : protocol information
    temperature | float : temperature [K]
    cpu | int : cpu to allocate
    log : logger instance

    return None
    """

    log.info(f"STARTING PROTOCOL {p.number}")
    log.info(
        f"\nActive conformers for this phase: {len([i for i in conformers if i.active])}\n"
    )
    count = 1
    for i in conformers:
        if not i.active:
            continue
        if i.energies.get(str(p.number)):
            continue
        launch(count, i, p, cpu, log, temperature, conformers)
        count += 1

    conformers = sorted(conformers)

    calculate_rel_energies(conformers, temperature)

    log.info("\nEnded Calculations\n")
    create_summary("Summary", conformers, log)
    log.info(
        "\nTotal elapsed time: "
        + str(
            datetime.timedelta(
                seconds=sum(
                    [i._last_energy["time"] for i in conformers if i.active]
                )
            )
        )
    )

    log.debug("Start Pruning")
    conformers = check_ensemble(conformers, p, log)
    save_snapshot(f"ensemble_after_{p.number}.xyz", conformers, log)

    create_summary("Summary After Pruning", conformers, log)

    log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')

    return None


def last_protocol_completed(conf, idx: int) -> bool:
    """
    Getting if all the conformers have been calculated until the idx-th protocol

    conf | list : whole ensemble list
    idx | int : index of the last protocol executed

    return | bool
    """

    tmp = []
    for i in conf:
        if i.energies.get(int(idx)) is not None and i.active:
            tmp.append(i)

    return (
        len(
            [
                # i for i in conf if i.energies.get(int(idx)) is not None and i.active
                tmp
            ]
        )
        == 0
    )


def create_summary(title, conformers, log):
    """
    Create the summary of a ensemble information

    title | str : title of the summary
    conformers | list : whole ensemble list
    log : logger instance
    """

    log.info(title)
    log.info("")
    log.info(
        tabulate(
            [i.create_log() for i in conformers if i.active],
            headers=[
                "Conformers",
                "E[Eh]",
                "G[Eh]",
                "B[cm-1]",
                "E. Rel [kcal/mol]",
                "Pop [%]",
                "Elap. time [sec]",
            ],
            floatfmt=".6f",
        )
    )
    log.info("")

    return None


def start_calculation(
    conformers, protocol, cpu: int, temperature: float, start_from: int, log
) -> None:
    """
    Main calculation loop

    conformers | list : whole ensemble list
    protocol | list : whole protocol steps
    cpu | int : cpu allocated
    temperature | float : temperature [K]
    start_from | int : index of the last protocol executed
    log : logger instance

    return None
    """
    if start_from != 0:
        if last_protocol_completed(conformers, start_from):
            conformers = check_ensemble(conformers, protocol[start_from], log)
            create_summary("Summary", conformers, log)

    for p in protocol[start_from:]:
        with open("last_protocol", "w") as f:
            f.write(str(p.number))
        run_protocol(conformers, p, temperature, cpu, log)
        if p.graph:
            Graph(conformers, p, log, temperature)

    save_snapshot("final_ensemble.xyz", conformers, log)
    log.info(f'{"="*15}\nCALCULATIONS ENDED\n{"="*15}\n\n')

    calculate_rel_energies(conformers, temperature)
    create_summary("Final Summary", conformers, log)

    return None


def restart() -> tuple:
    """
    Reload ensemble, protocol and setting from previous calculation

    return | list, list, int
    """
    confs = json.load(open("checkpoint.json"))
    ensemble = [Conformer.load_raw(confs[i]) for i in confs]

    p = json.load(open("protocol_dump.json"))
    protocol = [Protocol(**p[i]) for i in p]

    with open("last_protocol") as f:
        start_from = int(f.readlines()[0])

    return ensemble, protocol, start_from


def create_protocol(p, log) -> list:
    """
    Create the steps for the protocol to be executed

    p | dict : JSON read file of the protocol
    thrs | dict : JSON read file of the thresholds
    log : logger instance

    return | list : protocol steps
    """

    protocol = []

    log.info("Loading Protocol\n")
    last_prot_with_freq = None
    for idx, d in p.items():
        func = d.get("functional", None)
        add_input = d.get("add_input", "")
        graph = d.get("graph", False)

        if not func:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nFUNCTIONAL key must be passed in order to calculate energy.\nDFT functional or HF for Hartree-Fock calculation or semi-empirical methods (XTB1/XTB2/PM3/AM1 or similar supported by the calculator) (Problem at {ordinal(int(idx))} protocol definition)\n{'='*20}\nExiting\n{'='*20}\n"
            )
            raise IOError(
                "There is an error in the input file with the definition of the functional. See the output file."
            )

        if graph:
            if not add_input:
                log.critical(
                    f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nADD_INPUT must be set so the proper calculation (TD-DFT or CASSCF/RASSCF) to simulate the electronic spectra (Problem at {ordinal(int(idx))} protocol definition)\n{'='*20}\n Exiting"
                )
                raise IOError(
                    "There is an error in the input file with the definition of the functional. See the output file."
                )
            if not last_prot_with_freq:
                log.critical(
                    f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nElectrical spectra requires Boltzmann population over ∆G. In the specified protocol there is NO frequency calculation turned on.\n{'='*20}\nExiting\n{'='*20}"
                )
            else:
                # check_protocol_grapher()
                pass

        if not graph and d.get("freq"):
            last_prot_with_freq = int(idx)

        protocol.append(Protocol(number=idx, **d))

    log.info(
        "\n".join(
            (
                f"{i.number}: {str(i)} - {i.calculation_level}\n {i.thr}"
                for i in protocol
            )
        )
        + "\n"
    )
    return protocol


def main():
    """
    Main script loop
    """

    args = parser_arguments()

    # Trying to reload the damped setting from a previously calculation. Else damp the settings

    if os.path.exists(os.path.join(os.getcwd(), "settings.json")):
        settings = json.load(open("settings.json"))
    else:
        settings = {
            "output": args.output,
            "cpu": args.cpu,
            "temperature": args.temperature,
        }
        json.dump(settings, open("settings.json", "w"), indent=4)

    # create the setting dictionary
    output = (
        settings.get("output", args.output)
        if not args.restart
        else ".".join(settings.get("output", args.output).split(".")[:-1])
        + "_restart.out"
    )
    cpu = settings.get("cpu", args.cpu)
    temperature = settings.get("temperature", args.temperature)

    # initiate the log
    log = create_log(output)

    if args.restart:
        # reload the previous information from checkpoint file
        conformers, protocol, start_from = restart()

    else:
        protocol = create_protocol(load_protocol(args.protocol), log)
        start_from = protocol[0].number
        json.dump(
            {i.number: i.__dict__ for i in protocol},
            open("protocol_dump.json", "w"),
            indent=4,
            cls=SerialiseEncoder,
        )
        conformers = read_ensemble(
            args.ensemble, args.charge, args.multiplicity, log
        )

    # start the loop
    start_calculation(
        conformers=conformers,
        protocol=protocol,
        cpu=cpu,
        temperature=temperature,
        start_from=int(start_from),
        log=log,
    )
