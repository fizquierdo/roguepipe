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
# cp testdata/MockNaRok/* testdata/   # uncomment this line to test non corvengence

# get the two best trees (constrained, normal)

echo "trying to find problematic taxa "


iterFolder=iteration_no$(($iter_id-1)) 

normMLTree=$( ls ./$iterFolder/search_normal/RAxML_bestTree*)
constMLTree=$( ls ./$iterFolder/search_normal/RAxML_bestTree*)

if [  $(ls $normMLTree | wc -l  ) -ne 1 -o $(ls $constMLTree | wc -l  ) -ne 1 ] ; then 
    echo "could not find best trees from constraint and unconstraint searches "
    exit 
fi 

cat $normMLTree $constMLTree > bothTrees

./scripts/newConstraintMaster.sh $cons_tree  bothTrees 

# overwrite old stuff 
\cp new_constraint/problematic.txt blacklist.txt
\cp new_constraint/newConstraint.txt  $cons_tree



if [ $(cat $cons_tree | wc -l  ) -ne 1 ] ; then
    echo "could not create new constraint tree "
    exit
fi

blacklist=blacklist.txt

# some cleanup 
rm -rf new_constraint
rm bothTrees

if [ $(cat $blacklist | wc -l ) -eq 0  ]; then
    echo "could NOT find any incongruent taxa. Stopping the search."
    exit
fi


# =======================================================

# Overwrite phylip file in $aln according to the blacklist from Andre rogues
ruby scripts/phylip_pruner.rb $blacklist $aln $aln
rm $blacklist

# Call myself again with pruned dataset and constraint
$0 $iter_id $aln $cons_tree 

