#!/bin/sh

# Remove results from other runs
rm -rf iteration_no*
rm CONSEL_ITERATIONS
# Restore original testdata
rm testdata/blacklist.txt
cp testdata/original/* testdata/

# Start from iteration 1 
./treepruner.sh 1 testdata/aln.phy testdata/cons_tree.nw 
