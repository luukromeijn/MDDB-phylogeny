from Bio import Phylo
from Bio import SeqIO

ref_tree = Phylo.read('data\RAxML_bestTree.output_1.tre', 'newick')
test_tree = Phylo.read('data\subtree.tre', 'newick')
new_clade = test_tree.root
Phylo.draw(ref_tree)
for clade in ref_tree.find_clades():
    if clade.name == 'Fi1cons_02':
        clade.name = new_clade.name
        clade.clades = new_clade.clades
        clade.branch_length = new_clade.branch_length   
Phylo.draw(ref_tree)

# data = SeqIO.index('Reference_ITS.fasta', 'fasta')
# for line in data:
#     print(line, data[line].seq[:10] + '...' + data[line].seq[-10:])
#     # print(data[line].seq)
