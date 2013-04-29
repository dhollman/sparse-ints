#!/usr/bin/env zsh
rsync -rPh --no-motd --inplace -u *.cc *.h *.py hopper:sparse_ints/

