code used to generate figures for publication

## how to run

first run `analysis.r` to generate initial analysis outputs (which get saved to the `out` subdirectory).
then run `figures.r` to generate figures, which generates a list of all relevant figures in a variable named `fig`.

assumes presence of cellranger outputs and associated metadata used for experiments under the `data` subdirectory.

## dependency management

`pak` was used for dependency management, and a lockfile is available (`pkg.lock`).