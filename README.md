# Gromos2Amber

Converts Gromos topology files to Amber topology files

## Usage

```
gromos2amber [-h]
             [--config_in INPUT_CONFIGURATION_FILE | --num_solvent N]
             [--config_out OUTPUT_CONFIGURATION_FILE]

optional arguments:
    -h, --help            show this help message and exit
    --config_in INPUT_CONFIGURATION_FILE
                          Input Gromos-format configuration file
    --num_solvent N       Number of solvent molecules to include in output
                          topology. (Default: 0)
    --config_out OUTPUT_CONFIGURATION_FILE
                          Output Amber-format configuration file
```

## Versions

First version number increments with significant changes that may impact usage.

Second version number increments with changes that might change the format
of the output but are not expected to impact simulation dynamics.

Third version number increments with any minor change that should not change
the output.
