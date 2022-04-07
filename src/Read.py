from Bio import Phylo
from Bio import SeqIO
from alfpy.utils import seqrecords, distmatrix
from alfpy.utils.data import subsmat
from alfpy import wmetric
from alfpy import ncd

def find_names(clade):
    return

ref_tree = Phylo.read('RAxML_bestTree.output_1.tre', 'newick')
# test_tree = Phylo.read('test.tre', 'newick')
# new_clade = test_tree.root
# Phylo.draw(ref_tree)
# for clade in ref_tree.find_clades():
#     if clade.name == 'Fi1cons_02':
#         clade.name = new_clade.name
#         clade.clades = new_clade.clades
#         clade.branch_length = new_clade.branch_length   
# Phylo.draw(ref_tree)
matrix = subsmat.get('blosum62')
unknown = open('Unknown_ITS.fasta')
unknown_records = seqrecords.read_fasta(unknown)
unknown.close()
fh = open('Reference_ITS.fasta')
seq_records = seqrecords.read_fasta(fh)
fh.close()
for x in range(1):
    seq_records.add(unknown_records.id_list[x], unknown_records.seq_list[x])
    # dist = ncd.Distance(seq_records)
    dist = wmetric.Distance(seq_records, matrix)
    distances = []
    for y in range(seq_records.count-1):
        distances.append((seq_records.id_list[y], dist.pairwise_distance(y, seq_records.count-1)))
    best_distances = sorted(distances, key = lambda i: i[1])[:5]
    # print(best_distances)
    ancestor = ref_tree.common_ancestor(best_distances[0][0], best_distances[1][0], best_distances[2][0], best_distances[3][0], best_distances[4][0])
    # find_names(ancestor)
    for clade in ancestor.clades:
        if clade.name:
            print(clade.name)
        else:
            print(clade.clades)
Phylo.draw(ref_tree)
    # for distance in best_distances:
    #     for clade in ref_tree.find_clades(name=distance[0]):
    #         print(clade.name)
    # for clade in ref_tree.find_clades(name=best_distances[0][0]):
    #     print(clade.clades)
# new_file = open('new_file.fasta', 'w')
# new_file.write(">" + seq_records.id_list[-1] + "\n" + seq_records.seq_list[-1] + "\n")
# for distance in best_distances:
#     sequence = seq_records.seq_list[seq_records.id_list.index(distance[0])]
#     new_file.write(">" + distance[0] + "\n" + sequence + "\n")
# new_file.close()

# data = SeqIO.index('Reference_ITS.fasta', 'fasta')
# for line in data:
#     print(line, data[line].seq[:10] + '...' + data[line].seq[-10:])
#     # print(data[line].seq)
