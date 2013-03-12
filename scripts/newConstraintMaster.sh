#! /bin/bash


if [ $# != 2 ]; then
    echo "./script currentConstraintTree bothTrees "
    exit
fi

mkdir -p  new_constraint
cd  new_constraint


../scripts/getProblematicTaxa.sh ../$2 loose  > problematic.txt 

if [ "$(cat problematic.txt | wc -l )" == "0" ] ; then
    rm problematic.txt 
    ../scripts/getProblematicTaxa.sh ../$2 thorough  > problematic.txt 
fi 

cat ../$1 | tr "," "\n" | tr -d  "();"  >  currentTaxa

../scripts/getNewConstraintTree.py currentTaxa problematic.txt | grep -v module  2>/dev/null  > newConstraint.txt 

cd .. 
