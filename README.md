# Conformer Ensemble Pruning Software

This software is a tool designed for molecular modeling to simplify and optimize their pruning and optimizing workflow of an ensemble of conformers. 

With this software, users can quickly and easily identify the most important conformers and reduce the computational burden of simulations.

---

## Installation

To use the software, follow these steps:

- Clone the repository onto your local machine using Git.
- Install the required dependencies listed in the requirements.txt file. This can be done by running the following command in your terminal:
```bash
pip install -r requirements.txt
```

## Usage

To use the software, first, create two JSON files specifying the protocol to use and conformers' pruning parameters. 
The file format is described in detail
```bash
python ensemble_analyser.py -h-p # to obtain an example of the protocol.json file
python ensemble_analyser.py -h-t # to obtain an example of the treshold.json file
```

Next, run the `ensemble_analyser.py` script with the path to the JSON file as a command-line argument. The script will read the input file, perform the conformer pruning, and output a file with the pruned conformers.

```bash
python ensemble_analyser.py ensemble.xyz
````

## Parameters

The software uses the following parameters to determine the most relevant conformers:

- thrG: Refers to the energy (G or E) threshold to consider two conformers equivalent together with thrB.
- thrB: Refers to the rotary constant threshold to consider two conformers equivalent together with thrG.
- thrGMAX: Refers to the maximum energy window considered. Conformers lying above it will be sorted out immediately.

---

## Contributing

Contributions to this software are always welcome! If you have any ideas or suggestions, please feel free to submit a pull request or open an issue.

## License

This software is licensed under the GPL-3.0 License. See the LICENSE file for details.
