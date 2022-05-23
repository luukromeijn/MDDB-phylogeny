from Supertree import Backbone, RepresentativesTree

# Generating the tree #TODO saving the entire file still doesn't work
rep_tree = RepresentativesTree('results/supertree/representatives_tree.tre')
backbone = Backbone(representatives_tree=rep_tree)

# Displaying the tree (local device only)
# backbone = Backbone('src/backbone_with_taxon.tre')
# backbone.show()
