from Bio import Phylo
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from alfpy.utils import seqrecords
from alfpy.utils.data import subsmat
from alfpy import wmetric
from alfpy import word_pattern
from alfpy import word_vector
from alfpy import word_distance
import os
import glob
from tqdm import tqdm

test_tree = Phylo.read('RAxML_bestTree.test_tree.tre', 'newick')
Phylo.draw(test_tree)
kmer_ref_tree = Phylo.read('kmer_ref_tree.newick', 'newick')
Phylo.draw(kmer_ref_tree)
kmer_raxml_tree = Phylo.read('kmer_RAxML_tree.newick', 'newick')
Phylo.draw(kmer_raxml_tree)
wmetric_ref_tree = Phylo.read('wmetric_ref_tree.newick', 'newick')
Phylo.draw(wmetric_ref_tree)
wmetric_raxml_tree = Phylo.read('wmetric_RAxML_tree.newick', 'newick')
Phylo.draw(wmetric_raxml_tree)

def find_names(ancestor, names):
    for clade in ancestor.clades:
        if clade.name:
            names.append(clade.name)
        else:
            find_names(clade, names)

def write_to_fasta(childs, records):
    new_file = open('new_file.fasta', 'w')
    new_file.write(">" + records.id_list[-1] + "\n" + records.seq_list[-1] + "\n")
    for child in childs:
        sequence = records.seq_list[records.id_list.index(child)]
        new_file.write(">" + child + "\n" + sequence + "\n")
    new_file.close()

def update_tree():
    ref_tree = Phylo.read('RAxML_bestTree.output_1.tre', 'newick')
    # test_tree = Phylo.read('new_ref_tree.newick', 'newick')
    # Phylo.draw(ref_tree)
    # Phylo.draw(test_tree)

    # matrix = subsmat.get('blosum62')
    unknown = open('Unknown_ITS.fasta')
    unknown_records = seqrecords.read_fasta(unknown)
    unknown.close()
    fh = open('Reference_ITS.fasta')
    seq_records = seqrecords.read_fasta(fh)
    fh.close()
    for x in tqdm(range(unknown_records.count)):
        if unknown_records.seq_list[x] not in seq_records.seq_list:
            seq_records.add(unknown_records.id_list[x], unknown_records.seq_list[x])
            pattern = word_pattern.create(seq_records.seq_list, word_size=6)
            counts = word_vector.Counts(seq_records.length_list, pattern)
            dist = word_distance.Distance(counts, 'google')
            # dist = wmetric.Distance(seq_records, matrix)
            distances = []
            for y in range(seq_records.count-1):
                distances.append((seq_records.id_list[y], dist.pairwise_distance(y, seq_records.count-1)))
            best_distances = sorted(distances, key = lambda i: i[1])[:3]
            # print(best_distances)
            ancestor = ref_tree.common_ancestor(best_distances[0][0], best_distances[1][0], best_distances[2][0]) #, best_distances[3][0], best_distances[4][0]
            child_names = []
            find_names(ancestor, child_names)
            # print(len(child_names))
            write_to_fasta(child_names, seq_records)
            # Phylo.draw(ref_tree)

            # TRYING TO GET MAFFT AND RAXML TO WORK FROM WITHIN PYTHON BY CALLING COMMANDLINE
            mafft_exe = "mafft-win\mafft.bat"
            in_file = "new_file.fasta"
            mafft_cline = MafftCommandline(mafft_exe, input=in_file)
            stdout, stderr = mafft_cline()
            with open("aligned.fasta", "w") as handle:
                handle.write(stdout)
            raxml_exe = "raxmlHPC-PTHREADS-SSE3.exe"
            raxml_cline = RaxmlCommandline(raxml_exe, sequences="aligned.fasta", model="GTRGAMMA", name="new_tree")
            output, error = raxml_cline()
            new_subtree = Phylo.read('RAxML_bestTree.new_tree', 'newick')
            new_clade = new_subtree.root
            ancestor.clades = new_clade.clades
            # DELETE ALL RAXML FILES
            files = glob.glob('*new_tree')
            for f in files:
                os.remove(f)
        else:
            print('Zit er al in')
    Phylo.write(ref_tree, 'kmer_RAxML_tree.newick', 'newick')
    Phylo.draw(ref_tree)