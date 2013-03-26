#! /usr/bin/python

import sys
from taxonomy import * 



if(len(sys.argv) != 3 ) : 
    sys.stderr.write("USAGE: %s taxFile ncbi-taxonomy\n\n where taxFile is a file containing all relevant taxa (one per line)\n" %  sys.argv[0])
    sys.exit(0)

listOfTaxa= map(lambda x : x.strip(),open(sys.argv[1], "r").readlines())

tax = Taxonomy()
tax.init_parseTaxFile(sys.argv[2])
tax.reduceTaxonomy(listOfTaxa)
print tax.getNewickString() + ";"
