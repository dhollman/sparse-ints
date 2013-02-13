#!/usr/bin/env zsh
rsync -rPh --no-motd --inplace -u *.cc *.h hopper:sparse_ints/
rsync -rPh --no-motd --inplace -u /Users/dhollman/Projects/ao-mp2-r12/python/int_pics.py hopper:aoints/

