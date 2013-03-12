#! /usr/bin/python

from ete2 import * 
import sys 


class ExtendedTree(Tree):
    def getUnrootedLCA(self, listOfTaxa): 
        # not sure if this is posisble 
        assert( 8==0 )          # :-)
        
        oldRoot = tr 

        leaves = self.get_leaf_names()
        assert(set(listOfTaxa)  <= set(leaves))        
        outgroupTaxon = list(set(leaves) - set(listOfTaxa))[0]
        
        self.set_outgroup(outgroupTaxon)
        lca = self.get_common_ancestor(listOfTaxa)
        
        


    def getSubTreeLength(self,node):
        return node.dist +  sum(map(lambda x : x.dist, node.get_descendants()))

    def rootByInnermostNode(self): 
        part = self.getInnermostPartition()
        self.unroot()
        self.set_outgroup(self.get_common_ancestor(list(part)))

    
    def getInnermostPartition(self): 
        innerNodes = set(self.get_descendants()) -  set(self.get_leaves())
        treeLength = self.getSubTreeLength(self) 

        for n in self.get_descendants(): 
            n.add_feature("inScore", 0.0)

        bestNode = list(innerNodes)[0]
        bestScore = 0
        for n in innerNodes:             
            cs = n.get_children()
            l1 = self.getSubTreeLength(cs[0])
            l2 = self.getSubTreeLength(cs[1])

            self.set_outgroup(n) 
            toAdd = n.up.dist + n.dist
            cs = n.up.get_children()
            sis = n.get_sisters()
            assert(len(sis) == 1)
            assert(cs[0] == n or cs[1] == n)
            
            sis = sis[0]
            l3 = self.getSubTreeLength(sis)
            l3 += toAdd

            assert(treeLength - l1 - l2 - l3 < 0.00001 )
            score = treeLength - max(max(l1,l2),l3)

            if score > bestScore: 
                bestNode = n
                bestScore = score
        return bestNode.get_leaf_names()
            

if __name__ == "__main__": 
    if len(sys.argv) != 2: 
        print "need tree file"
        sys.exit()

    treeFile = sys.argv[1]
    tr = ExtendedTree(treeFile)
    # tr.getUnrootedLCA(["AlligatorMVZ"])
    print tr 
    
