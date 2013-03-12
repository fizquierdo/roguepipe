#! /bin/bash

if [ $# != 2 ]; then
    echo -e  "./script file [thorough|loose]\nwhere file is a newick file containing two trees and  thorough indicates, how aggressively taxa shall be removed (try loose first). "
    exit 
fi

relPath=$(dirname $0)

file=$1

runid=$RANDOM

mkdir -p rnr_run
cd rnr_run
rnr="$(pwd)/../../bin/RogueNaRok"

$rnr -i ../$file -c 100 -s 2 -n $runid > /dev/null

tail -n +3  RogueNaRok_droppedRogues.$runid |  cut -f 3,4 > rogues 

while read line ; 
do 
    spec=$( echo $line  | cut -f 1 -d ' '  ) 
    harm=$( echo $line | cut -f 2 -d ' '  ) 

    if [ "$2" == "loose"  -a    "$(echo $harm | cut -f 1 -d '.')" -gt 1 ]; then
	echo "$spec" | tr "," "\n"
    elif [ "$2" == "thorough"  ] ; then 
	echo "$spec" | tr "," "\n"
    else
	break 
    fi
    
done  < rogues   

rm rogues 

cd .. 
