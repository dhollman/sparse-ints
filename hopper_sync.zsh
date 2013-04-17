#!/usr/bin/env zsh
rsync -rPh --no-motd --inplace -u *.cc *.h *.py edison.nersc.gov:sparse_ints/

