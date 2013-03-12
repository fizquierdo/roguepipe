#! /usr/bin/python

import sys 
from taxonomy import * 

if len(sys.argv) != 3 :
    print "./script oldListOfTaxa taxaToPrune" 
    sys.exit()

oldTaxa = set(map(lambda x : x.strip(),  open(sys.argv[1],"r").readlines())) 
removeTaxa = set(map(lambda x : x.strip(), open(sys.argv[2],"r").readlines())) 

taxaNowNeeded = list(oldTaxa - removeTaxa) 


assert(len(taxaNowNeeded) != 0)

tax = Taxonomy()

tax.init_parseTaxFile( "../ncbi.dmp.cleaned")

tax.reduceTaxonomy(taxaNowNeeded)

print tax.getNewickString()
