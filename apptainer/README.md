# README for PyPHLAWD, apptainer

- Last modified: 2025-11-03 11:49:07
- Sign: nylander

## Description

[PyPHLAWD.def](PyPHLAWD.def) - Apptainer definition file for
[PyPHLAWD](https://github.com/FePhyFoFum/PyPHLAWD/),
[phyx](https://github.com/FePhyFoFum/phyx), and
[phlawd_db_maker](https://github.com/blackrim/phlawd_db_maker).

## Build and run

Tested using apptainer v1.4.4 on Ubuntu 24.04.

    $ sudo apptainer build PyPHLAWD.sif PyPHLAWD.def

    $ ./PyPHLAWD.sif

## Notes

- Installing networkx using apt (`apt install python3-networkx`) installs
  an old version of networkx (see deprecation
  <https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations>)

