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

# kmer_raxml_tree = Phylo.read('Reference_tree.newick', 'newick')
# Phylo.draw(kmer_raxml_tree)
# # print("through raxml")
# test_tree = Phylo.read('20p_distance_OTU_tree.newick', 'newick')
# Phylo.draw(test_tree)
# print("kmer ref tree")
# kmer_ref_tree = Phylo.read('test_fi1cons_tree.newick', 'newick')
# Phylo.draw(kmer_ref_tree)
# # print("kmer raxml tree")
# Phylo.draw(kmer_raxml_tree)
# print("wmetric ref tree")
# wmetric_ref_tree = Phylo.read('wmetric_ref_tree.newick', 'newick')
# Phylo.draw(wmetric_ref_tree)
# # print("wmetric raxml tree")
# wmetric_raxml_tree = Phylo.read('wmetric_RAxML_tree.newick', 'newick')
# Phylo.draw(wmetric_raxml_tree)

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

def write_to_fasta_again(childs, records, extra):
    new_file = open('new_file.fasta', 'w')
    new_file.write(">" + extra + "\n" + records.seq_list[records.id_list.index(extra)] + "\n")
    for child in childs:
        sequence = records.seq_list[records.id_list.index(child)]
        new_file.write(">" + child + "\n" + sequence + "\n")
    new_file.close()

def create_subtree(ancestor):
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
    ancestor.name = new_clade.name
    # DELETE ALL RAXML FILES
    files = glob.glob('*new_tree')
    for f in files:
        os.remove(f)

def update_tree():
    # ref_tree = Phylo.read('Reference_tree.newick', 'newick')
    ref_tree = Phylo.read('RAxML_bestTree.partly_tree.tre', 'newick')
    # test_tree = Phylo.read('new_ref_tree.newick', 'newick')

    # matrix = subsmat.get('blosum62')
    unknown = open('Reference_ITS.fasta')
    unknown_records = seqrecords.read_fasta(unknown)
    unknown.close()
    fh = open('partly_fi1cons.fasta')
    seq_records = seqrecords.read_fasta(fh)
    fh.close()
    skipped = []
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
            sorted_distances = sorted(distances, key = lambda i: i[1])
            best_distances = [sorted_distances[x][0] for x in range(len(sorted_distances)) if sorted_distances[x][1] <= 1.1*sorted_distances[0][1]]
            print(unknown_records.id_list[x], best_distances)
            if len(best_distances) > 2:
                extra_seqs = []
                for skipped_seq in skipped:
                    if skipped_seq in best_distances:
                        extra_seqs.append(skipped_seq)
                        best_distances.remove(skipped_seq)
                for extra in extra_seqs:
                    skipped.remove(extra)
                ancestor = ref_tree.common_ancestor(best_distances)
                child_names = []
                if len(ancestor.clades) == 0:
                    child_names.append(ancestor.name)
                else:
                    find_names(ancestor, child_names)
                child_names.extend(extra_seqs)
                write_to_fasta(child_names, seq_records)
                create_subtree(ancestor)
                # Phylo.draw(ref_tree)
            else:
                skipped.append(unknown_records.id_list[x])
        else:
            print('Zit er al in')
    if skipped:
        print(len(skipped))
        for skip in skipped:
            names = [terminal.name for terminal in ref_tree.get_terminals()]
            if skip not in names:
                new_distances = []
                for a in range(seq_records.count):
                    new_distances.append((seq_records.id_list[a], dist.pairwise_distance(a, seq_records.id_list.index(skip))))
                new_sorted_distances = sorted(new_distances, key = lambda i: i[1])[1:4]
                best_new_dist = [new_sorted_distances[x][0] for x in range(len(new_sorted_distances))]
                extra_seqs = []
                for skipped_seq in skipped:
                    if skipped_seq in best_new_dist:
                        extra_seqs.append(skipped_seq)
                        best_new_dist.remove(skipped_seq)
                ancestor = ref_tree.common_ancestor(best_new_dist)
                child_names = []
                find_names(ancestor, child_names)
                child_names.extend(extra_seqs)
                write_to_fasta_again(child_names, seq_records, skip)
                #BUILD SUBTREE
                create_subtree(ancestor)
    Phylo.write(ref_tree, 'test_fi1cons_tree.newick', 'newick')
    Phylo.draw(ref_tree)

# update_tree()