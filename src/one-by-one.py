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
from alfpy.utils import distmatrix
import os
import glob
from tqdm import tqdm

# kmer_raxml_tree = Phylo.read('single_alignment.tre', 'newick')
# Phylo.draw(kmer_raxml_tree)
# print("through raxml")
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
    for child in childs:
        sequence = records.seq_list[records.id_list.index(child)]
        new_file.write(">" + child + "\n" + sequence + "\n")
    new_file.write(">" + records.id_list[-1] + "\n" + records.seq_list[-1] + "\n")
    new_file.close()

def write_to_fasta_again(childs, records, extra):
    new_file = open('new_file.fasta', 'w')
    for child in childs:
        sequence = records.seq_list[records.id_list.index(child)]
        new_file.write(">" + child + "\n" + sequence + "\n")
    new_file.write(">" + extra + "\n" + records.seq_list[records.id_list.index(extra)] + "\n")
    new_file.close()

def create_subtree(ancestor, outgroups):
    # TRYING TO GET MAFFT AND RAXML TO WORK FROM WITHIN PYTHON BY CALLING COMMANDLINE
    mafft_exe = "mafft-win\mafft.bat"
    in_file = "new_file.fasta"
    mafft_cline = MafftCommandline(mafft_exe, input=in_file, inputorder = False) #, adjustdirection=True
    stdout, stderr = mafft_cline()
    with open("aligned.fasta", "w") as handle:
        handle.write(stdout)
    raxml_exe = "raxmlHPC-PTHREADS-SSE3.exe"
    constraint = "constraint.newick"
    outgroups_string = ','.join(outgroups)
    if len(outgroups):
        raxml_cline = RaxmlCommandline(raxml_exe, sequences="aligned.fasta", binary_constraint=constraint, outgroup=outgroups_string, model="GTRCAT", name="new_tree")#, binary_constraint=constraint
    else:
        raxml_cline = RaxmlCommandline(raxml_exe, sequences="aligned.fasta", binary_constraint=constraint, model="GTRCAT", name="new_tree")#, binary_constraint=constraint
    output, error = raxml_cline()
    new_subtree = Phylo.read('RAxML_bestTree.new_tree', 'newick')
    # new_ancestor = new_subtree.common_ancestor(outgroups)
    for outgroup in outgroups:
        new_subtree.prune(outgroup)
    new_clade = new_subtree.root
    ancestor.clades = new_clade.clades
    ancestor.name = new_clade.name
    # DELETE ALL RAXML FILES
    files = glob.glob('*new_tree')
    for f in files:
        os.remove(f)

def reversal():
    fasta_sequences = SeqIO.parse(open('Unknown_ITS.fasta'),'fasta')
    with open('new_unknown.fasta', 'w') as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            out_file.write('>' + name + '\n')
            new_sequence = sequence[::-1]
            out_file.write(new_sequence + '\n')

def get_outgroups(ancestor, tree):
    path = tree.get_path(ancestor)
    if len(path):
        outgroup_parent = tree.from_clade(path[-3])
        ancestor_childs = [terminal.name for terminal in ancestor.get_terminals()]
        terminals = [terminal.name for terminal in outgroup_parent.get_terminals() if terminal.name not in ancestor_childs]
        for target in terminals[2:]:
            outgroup_parent.prune(target)
        return terminals[:2], outgroup_parent
    else:
        return [], tree

def update_tree():
    ref_tree = Phylo.read('single_alignment.tre', 'newick')
    # ref_tree = Phylo.read('RAxML_bestTree.partly_tree.tre', 'newick')
    # test_tree = Phylo.read('Reference_tree.newick', 'newick')

    # unknown = open('new_unknown.fasta')
    unknown = open('876_Ascomycota_Leotiomycetes_Helotiales.fasta')
    unknown_records = seqrecords.read_fasta(unknown)
    unknown.close()
    fh = open('all_in_one_unaligned.fasta')
    seq_records = seqrecords.read_fasta(fh)
    fh.close()
    count = 0
    avg_amount = 0
    avg_childs = 0
    for x in tqdm(range(unknown_records.count)):
        if unknown_records.seq_list[x] not in seq_records.seq_list:
            count += 1
            distances = []
            seq_records.add(unknown_records.id_list[x], unknown_records.seq_list[x])
            pattern = word_pattern.create(seq_records.seq_list, word_size=6)
            counts = word_vector.Counts(seq_records.length_list, pattern)
            dist = word_distance.Distance(counts, 'google')
            for y in range(seq_records.count-1):
                distances.append((seq_records.id_list[y], dist.pairwise_distance(y, seq_records.count-1)))
            sorted_distances = sorted(distances, key = lambda i: i[1])
            best_distances = [sorted_distances[x][0] for x in range(len(sorted_distances)) if sorted_distances[x][1] <= 1.025*sorted_distances[0][1]]
            if len(best_distances):
                avg_amount += len(best_distances)
                ancestor = ref_tree.common_ancestor(best_distances)
                outgroups, constraint = get_outgroups(ancestor, ref_tree)
                Phylo.write(constraint, 'constraint.newick', 'newick')
                child_names = []
                if len(ancestor.clades) == 0:
                    child_names.append(ancestor.name)
                else:
                    find_names(ancestor, child_names)
                child_names.extend(outgroups)
                avg_childs += len(child_names)
                write_to_fasta(child_names, seq_records)
                create_subtree(ancestor, outgroups)
        else:
            print('Zit er al in')
            # pass
    print(avg_amount/count)
    print(avg_childs/count)
    Phylo.write(ref_tree, 'backbone_plus_Ascomycota_Leotiomycetes_Helotiales.newick', 'newick')
    Phylo.draw(ref_tree)

update_tree()
# create_subtree()
# reversal()