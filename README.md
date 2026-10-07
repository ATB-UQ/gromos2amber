# gromos2amber

Converts GROMOS topology (and optionally coordinate) files into AMBER `prmtop` / `inpcrd` files.
Pure Python, no third-party dependencies. **Status: live-support** -- used by the website's system
topology download and by the AMBER/pmemd paths of the condensed-phase and benchmark work.

## Usage

Python API (what the platform uses):

```python
from gromos2amber import convert, GromosFormatError, IllegalArgumentError
convert("topology.top", "topology.prmtop",
        config_in="positions.g96", config_out="positions.inpcrd",
        solvent_resname="SOL", num_solvent=-1)
```

Command line, via the script `bin/gromos2amber` (**not installed by pip**: `pyproject.toml`
declares no `scripts`/entry point, so `.venv/bin/gromos2amber` does not exist; run
`python bin/gromos2amber ...` or add an entry point):

```
gromos2amber [--config_in G96 | --num_solvent N] [--config_out INPCRD]
             [--solvent_resname SOL] < topology.top > topology.prmtop
```

`--config_in` and `--num_solvent` are mutually exclusive: with a configuration the solvent
count is taken from it, otherwise `--num_solvent` solvent molecules are included (default 0).

## Layout

`gromos2amber/`: `Converter.convert` (entry), `GromosTopologyParser`, `Topology`, `Configuration`,
`AmberTopologyWriter`, `AmberConfigurationWriter`, `fortran_format`, `gromos_format`, `Errors`.
Note GROMOS quartic bond/angle terms are refitted to harmonic on output.

## Platform fit

- Used by: `website` (`site_pages/system_topology_download.py`), `gromos_job_wrapper`
  (`gromos2amber` job option, `deployment/gadi_pmemd.sh` / Setonix pmemd shim),
  `atb_condensed_phase` (`engines/amber.py`), `ncmas_benchmarks` (`scripts/build_inputs.py`).
- Not related to `amber_converter` (an unfinished GROMACS/ParmEd prototype; see its README).
- `siblings.txt` tracks branch `pip-installable-packaging`; version in `VERSION` (1.2.1).

## Tests

`tests/` is empty in this checkout (a tracked gitlink with no content), so there are no tests.
Validate a change by converting a molecule and loading the result with AMBER/ParmEd.
