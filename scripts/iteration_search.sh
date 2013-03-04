#!/bin/bash
if [ $# -ne 5 ]
then
  echo "Usage: `basename $0` basedir name aln tree_constain result_output_file"
  exit
fi
basedir=$1
name=$2
aln=$basedir/$3
cons_tree=$basedir/$4
result_name=$5
echo "Starting iteration $name with alignment ${aln} and constrain ${cons_tree}"

# Config raxml
num_threads=4
seed=12345
num_searches=3
rax=raxmlHPC-PTHREADS-SSE3
# Config consel
catpv=$basedir/conselbin/catpv
consel=$basedir/conselbin/consel
makermt=$basedir/conselbin/makermt

# Generate new folders
dir_iteration=$basedir/$name
dir_search_normal=${dir_iteration}/search_normal
dir_search_constrained=${dir_iteration}/search_constrained
dir_consel=${dir_iteration}/consel
mkdir $dir_iteration
mkdir $dir_search_normal
mkdir $dir_search_constrained
mkdir $dir_consel

# Constrained search
cd $dir_search_constrained
$rax -T $num_threads -n search_cons_${name} -s $aln -p $seed -m GTRGAMMA -g $cons_tree -N $num_searches 
cat RAxML_result.* > $dir_consel/${num_searches}_bunch_constrained.nw

# Normal search
cd $dir_search_normal
$rax -T $num_threads -n search_norm_${name} -s $aln -p $seed -m GTRGAMMA -N $num_searches 
cat RAxML_result.* > $dir_consel/${num_searches}_bunch_normal.nw

cd $dir_consel
# per-site LH
bunch=bunch_norm_cons.nw
cat ${num_searches}_bunch_normal.nw ${num_searches}_bunch_constrained.nw > $bunch
$rax -T $num_threads -f g -z $bunch -n $name -s $aln -m GTRGAMMA 
# Consel analysis
log=consel.log
cp RAxML_perSiteLLs.${name} ${name}.siteLH
echo "MAKERMT\n\n" >> $log
$makermt --puzzle ${name}.siteLH >> $log
echo "CONSEL\n\n" >> $log
$consel ${name} >> $log
echo "Visualize results\n\n" >>$log
$catpv ${name} >> $log

# Now log must contain a table in the end like:
# rank item    obs     au     np |     bp     pp     kh     sh    wkh    wsh |
#    1    3    0.0  0.640  0.481 |  0.481  0.198      0  0.659      0  0.659 |
#    2    2    0.0  0.640  0.481 |  0.481  0.198      0  0.659      0  0.659 |
#    3    4    0.1  0.685  0.482 |  0.286  0.175  0.463  0.686  0.463  0.686 |
#    4    5    0.1  0.685  0.482 |  0.286  0.175  0.463  0.686  0.463  0.686 |
#    5    6    0.1  0.685  0.482 |  0.286  0.175  0.463  0.686  0.463  0.686 |
#    6    1    0.9  0.206  0.235 |  0.234  0.080  0.301  0.315  0.338  0.400 |


# from column item we can decide if we need a further iteration or not
# if the first 3 items all belong to [1,2,3] (normal search) then we need to continue, otherwise we can stop, because one of the constrained searches has been "mixed in" the good group.
ruby ${basedir}/scripts/consel_evaluator.rb $log $num_searches >> ${basedir}/${result_name}

