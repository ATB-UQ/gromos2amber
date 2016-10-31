# Gromos2Amber

Converts Gromos topology files to Amber topology files

## Usage

```
gromos2amber [-h]
             [--config_in INPUT_CONFIGURATION_FILE | --num_solvent N]
             [--config_out OUTPUT_CONFIGURATION_FILE]
             < INPUT_GROMOS_TOPOLOGY
             > OUTPUT_GROMOS_TOPOLOGY

optional arguments:
    -h, --help            show this help message and exit
    --config_in INPUT_CONFIGURATION_FILE
                          Input Gromos-format configuration file
    --num_solvent N       Number of solvent molecules to include in output
                          topology. (Default: 0)
    --config_out OUTPUT_CONFIGURATION_FILE
                          Output Amber-format configuration file
    --solvent_resname SOLVENT_RESIDUE_NAME
                          The name of the solvent residues. Maximum 4
                          characters. (Default: SOL)
```

