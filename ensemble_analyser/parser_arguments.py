import argparse, os
import sys, json



def print_help_protocol():

    example = json.dumps(
        {
            "0": {
                "func": "str: DEFINE THE DFT FUNCTIONAL",
                "basis": "str: DEFINE THE BASIS SET FOR THE CALCULATION. DEFAULT: def2-svp",
                "opt": "bool: TRUE IF WANT TO OPTIMISE. DEFAULT: False",
                "freq": "bool: TRUE IF WANT ANALYTICAL FREQUENCY CALCULATION. Defaul: False",
                "solv" : {
                    "solvent": 'str|null: NAME OF THE SOLVENT. IF GAS PHASE DEFINE AS NULL',
                    "smd" : 'bool: TRUE IF SMD MODEL EMPLOYED FOR IMPLICIT CALCULATION, ELSE CPCM USE',
                    '_comment' : 'SOLV KEYWORD CAN BE OMMITED AND A GAS PHASE CALCULATION WILL BE CARRIED.'
                    }, 
                "add_input" : "str: ADDITIONAL INPUT TO BE PASSED TO THE CALCULATOR. NO SANITY CHECK ON THIS BLOCK! USE WITH CARE. DEFAULT null"
                },
        }, indent=4
    ) 

    txt = f'The protocol file must be a JSON file and it can contains as many calculation parts as it is desidable. These calculations can be mere single point, optimisations, solo frequency or optimisation-frequency calculations. The JSON file is formatted as follow, note that ONLY "func" and "solvent" (if used the solv keyword) are mandatory:\n{example}\n'

    print(txt)
    sys.exit()


def print_help_threshold():
    with open(str(os.path.join(os.path.dirname(__file__), 'parameters_file','default_threshold.json')), 'r') as j:
        contents = json.loads(j.read())
    example = json.dumps(contents, indent=4)

    txt = f'With this JSON file the different pruning parameters can be change. USE WITH CARE!\n\nIn here all "type" of calculations possibile hade differents thrasholds:\n\t- thrG: refers to the energy (G or E) threshold to consider two conformers equivalntes togheter with thrB\n\t- thrB: refers to the rotary constant threshold to to consider two conformers equivalntes togheter with thrG\n\t- thrGMAX: refers to the max energy windows considered. Conformers lying above it will be sorted out immediately\n\nTwo conformers are equivalent if G_(CONFi)-G(CONFi-1)<thrG AND B_(CONFi)-B_(CONFi-1)<thrB.\n\nThe default parametes are:\n{example}\n\nIf you want to change only some of the parameters, ALL parameters must be passed.\nDO NOT CHANGE THE KEYWORDS'

    print(txt)
    sys.exit()






def parser_arguments():

    parser = argparse.ArgumentParser(add_help=False)

    # ensemble file
    input_group = parser.add_argument_group('Input Files')
    input_group.add_argument('ensemble', help='The ensemble file. Could be an xyz file (preferably) or other type parsable by OpenBabel')
    input_group.add_argument('-p', '--protocol', help='JSON file containg the computational protocol. Default: %(default)s', default=os.path.join(os.path.dirname(__file__), 'parameters_file','default_protocol.json'))
    input_group.add_argument('-t', '--threshold', help='JSON file containg the threshold divided by calculation type. Default: %(default)s', default=os.path.join(os.path.dirname(__file__), 'parameters_file','default_threshold.json'))

    system_group = parser.add_argument_group('System Parameters')
    system_group.add_argument('cpu', type=int, help='Define the number of CPU used by the calculations')
    system_group.add_argument('-c', '--calculator', help='Define the calculator to use. Default %(default)s', choices=['orca'], default='orca')


    other_group = parser.add_argument_group('Other Parameters')
    other_group.add_argument('-o', '--output', help='Define the output filename. Defaut: %(default)s', default='output.out')


    help_group = parser.add_argument_group('Get help')
    help_group.add_argument("-h", "--help", action="help", help="Show this help message")
    help_group.add_argument('-h-p', '--help-protocol', help='Get help to format correcly the protocol JSON-file', action='store_true')
    help_group.add_argument('-h-t', '--help-threshold', help='Get help to format correcly the threshold JSON-file', action='store_true')

    a = sys.argv
    if '-h-p' in a or '--help-protocol' in a:
        return print_help_protocol()
    
    if '-h-t' in a or '--help-threshold' in a:
        return print_help_threshold()

    return parser.parse_args()



if __name__ == '__main__':
    parser_arguments()
