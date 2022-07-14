from audioop import avg
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
import sys
import glob
from tqdm import tqdm

sys.setrecursionlimit(2000)
path = "l0.2_s3_4_1500_o2.0_a1_constr/l0.2_s3_4_1500_o2.0_a1_constr/chunks/unaligned"
dir_list = os.listdir(path)
dict_files = {}
for file_name in dir_list:
    open_file = open("l0.2_s3_4_1500_o2.0_a1_constr/l0.2_s3_4_1500_o2.0_a1_constr/chunks/unaligned/" + file_name)
    record = seqrecords.read_fasta(open_file)
    open_file.close()
    id = file_name.split("_")[0]
    dict_files[id] = record

def find_names(ancestor, names):
    for clade in ancestor.clades:
        if clade.name:
            names.append(clade.name)
        else:
            find_names(clade, names)

def write_to_fasta(childs, records, new_records, index):
    new_file = open('new_file.fasta', 'w')
    for child in childs:
        # try:
        sequence = records.seq_list[records.id_list.index(child)]
        # except:
        #     sequence = new_records.seq_list[new_records.id_list.index(child)]
        new_file.write(">" + child + "\n" + sequence + "\n")
    new_sequence = new_records.seq_list[index]
    new_file.write(">" + new_records.id_list[index] + "\n" + new_sequence + "\n")
    # new_file.write(">" + records.id_list[-1] + "\n" + records.seq_list[-1] + "\n")
    new_file.close()

def create_subtree(ancestor, outgroups):
    # TRYING TO GET MAFFT AND RAXML TO WORK FROM WITHIN PYTHON BY CALLING COMMANDLINE
    mafft_exe = "mafft-win\mafft.bat"
    in_file = "new_file.fasta"
    mafft_cline = MafftCommandline(mafft_exe, input=in_file, inputorder = False, localpair=True) #, adjustdirection=True
    stdout, stderr = mafft_cline()
    with open("aligned.fasta", "w") as handle:
        handle.write(stdout)
    raxml_exe = "raxmlHPC-PTHREADS-SSE3.exe"
    constraint = "constraint.newick"
    outgroups_string = ','.join(outgroups)
    if len(outgroups):
        raxml_cline = RaxmlCommandline(raxml_exe, sequences="aligned.fasta", binary_constraint=constraint, outgroup=outgroups_string, model="GTRCAT", name="new_tree")
    else:
        raxml_cline = RaxmlCommandline(raxml_exe, sequences="aligned.fasta", binary_constraint=constraint, model="GTRCAT", name="new_tree")
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

def get_outgroups(ancestor, tree):
    path = tree.get_path(ancestor)
    # print(tree.root)
    if len(path):
        if ancestor.name:
            outgroup_parent = tree.from_clade(path[-3])
        elif len(path) == 1:
            outgroup_parent = tree.root
        else:
            outgroup_parent = tree.from_clade(path[-2])
        # Phylo.draw(test)
        ancestor_childs = [terminal.name for terminal in ancestor.get_terminals()]
        terminals = [terminal.name for terminal in outgroup_parent.get_terminals() if terminal.name not in ancestor_childs]
        for target in terminals[2:]:
            outgroup_parent.prune(target)
        return terminals[:2], outgroup_parent
    else:
        return [], tree

def update_tree():
    ref_tree = Phylo.read('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/supertree/backbone.tre', 'newick')
    # ref_tree = Phylo.read('single_alignment.tre', 'newick')

    # unknown = open('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/discarded/unidentified.fasta')
    # unknown = open('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/discarded/smallchunks.fasta')
    unknown = open('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/discarded/distant.fasta')
    # unknown = open('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/discarded/short_or_long.fasta')
    # unknown = open('mycodiversity_9_6_2022_1221.fa')
    # unknown = open('226_Ascomycota_Sordariomycetes_Sordariales.fasta')
    unknown_records = seqrecords.read_fasta(unknown)
    unknown.close()
    fh = open('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/supertree/representatives_unaligned.fasta')
    # fh = open('all_in_one_unaligned.fasta')
    seq_records = seqrecords.read_fasta(fh)
    fh.close()
    overall = open('l0.2_s3_4_1500_o2.0_a0_constr/l0.2_s3_4_1500_o2.0_a0_constr/supertree/backbone.fasta')
    overall_records = seqrecords.read_fasta(overall)
    overall.close()
    anc_dict = {}
    avg_childs = 0
    avg_dist = 0
    avg_smallest_dist = 0
    avg_worst_dist = 0
    avg_best = 0
    # for x in tqdm(range(unknown_records.count)):
    for x in tqdm(range(100)):
        if unknown_records.seq_list[x] not in seq_records.seq_list:
            distances = []
            seq_records.add(unknown_records.id_list[x], unknown_records.seq_list[x])
            pattern = word_pattern.create(seq_records.seq_list, word_size=6)
            counts = word_vector.Counts(seq_records.length_list, pattern)
            dist = word_distance.Distance(counts, 'google')
            for y in range(seq_records.count-2):
                distances.append((seq_records.id_list[y], dist.pairwise_distance(y, seq_records.count-1)))
            seq_records.id_list.remove(unknown_records.id_list[x])
            seq_records.seq_list.remove(unknown_records.seq_list[x])
            seq_records.count -= 1
            sorted_distances = sorted(distances, key = lambda i: i[1])
            best_distances = [sorted_distances[x][0] for x in range(len(sorted_distances)) if sorted_distances[x][1] <= 1.025*sorted_distances[0][1]]
            chunks = set()
            chunk_distances = []
            for name in best_distances:
                chunk = name.split("_SH")[0]
                chunks.add(chunk)
            for code in chunks:
                records = dict_files[code]
                records.add(unknown_records.id_list[x], unknown_records.seq_list[x])
                pattern = word_pattern.create(records.seq_list, word_size=6)
                counts = word_vector.Counts(records.length_list, pattern)
                dist = word_distance.Distance(counts, 'google')
                for y in range(records.count-2):
                    chunk_distances.append((records.id_list[y], dist.pairwise_distance(y, records.count-1)))
                records.id_list.remove(unknown_records.id_list[x])
                records.seq_list.remove(unknown_records.seq_list[x])
                records.count -= 1
            sorted_distances = sorted(chunk_distances, key = lambda i: i[1])
            best_distances = [sorted_distances[x][0] for x in range(len(sorted_distances)) if sorted_distances[x][1] <= 1.025*sorted_distances[0][1]]
            if len(best_distances):
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
                if len(child_names) < 500:
                    write_to_fasta(child_names, overall_records, unknown_records, x)
                    create_subtree(ancestor, outgroups)
                    overall_records.add(unknown_records.id_list[x], unknown_records.seq_list[x])

update_tree()