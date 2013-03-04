#!/bin/bash

if [ $# -ne 3 ]
then
  echo "Usage: `basename $0` iter_id aln tree_constain"
  exit
fi
iter_id=$1
aln=$2
cons_tree=$3

echo "Starting iteration $name with alignment ${aln} and constrain ${cons_tree}"

# Search iteration
result_name=CONSEL_ITERATIONS
scripts/iteration_search.sh $PWD iteration_no${iter_id} $aln $cons_tree $result_name
iter_result=`tail -n 1 $result_name`

# Stop or Prune and launch next iteration
if [ $iter_result = SIMILAR_TOPOLOGIES ]
then
  echo "DONE"
  echo "Constrained and unconstrained topologies are not significantly different after iteration ${iter_id}"
  exit
fi

# Next iteration
let iter_id+=1
echo "next iteration: ${iter_id}"
# Possibly we want to force a maximum number of iterations 
# (a manual restart is straight-forward)
max_iter=3
if [ $iter_id -gt $max_iter ]
then
  echo "Too many iterations. Aborting."
  exit
fi

# =======================================================
# TODOANDRE 
# Call rogue narok, generate rogue list and overwrite cons_tree
# $cons_tree, testdata/blacklist.txt <- narok($aln, $cons_tree)  

#  Mock Andres magic (change next line by real implementation)
cp testdata/MockNaRok/* testdata/   # uncomment this line to test non corvengence
# =======================================================

# Overwrite phylip file in $aln according to the blacklist from Andre rogues
ruby scripts/phylip_pruner.rb testdata/blacklist.txt $aln $aln

# Call myself again with pruned dataset and constraint
$0 $iter_id $aln $cons_tree 

