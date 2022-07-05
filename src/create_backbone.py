'''Creates a fungal phylogenetic backbone with data from UNITE.'''

import os, glob, sys
import numpy as np
import time
from Bio import Phylo
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Align.Applications import MafftCommandline
from Chunks import Chunk, chunk_size_report, discard_small_chunks, UniteData
from DistanceData import PairwiseDistance
from Supertree import Backbone, RepresentativesTree

def divide_data(result_dir: str, fasta_path: str, taxon_path: str, matrix_file_path: str,
    mafft_path: str, raxml_path: str, length_tolerance: float, min_split_depth: int, 
    max_split_depth: int, max_chunk_size: int, outlier_strictness: float, representatives_alg: int, 
    full_constraint: bool, localpair: bool):
    """Divides a UNITE dataset into smaller chunks.
    Builds an overarching tree containing 2 representatives for each chunk. 
    Identifies outgroups for all chunks.

    Parameters
    ----------
    result_dir : str
        Path to desired output directory
    fasta_path : str
        Path to FASTA file with UNITE data
    taxon_path: str
        Path to corresponding UNITE taxonomy file
    matrix_file_path: str
        Path to distance matrix file
    mafft_path: str
        Path to mafft executable
    raxml_path: str
        Path to raxml executable
    length_tolerance: float
        Allowed factor of deviation compared to the average sequence length in UNITE.
    min_split_depth: int
        The data is at least split up at this taxonomy rank (1=phylum, ... , 6=species)
    max_split_depth: int
        The data will not be splitted beyond this taxonomy rank.
    max_chunk_size: int
        Size threshold for which the algorithm will split chunks further (if allowed by max_split_depth).
    outlier_strictness: float
        Number of standard deviations that a sequences average distance to its chunk must be lower than 
        the average pairwise distance in the entire dataset.
    representatives_alg: int
        Determines which function is used to determine the chunk representatives. 
        If 0: pick sequences most average to their chunks.
        If 1: pick sequences most distant to each other within their chunk. 
    full_constraint: bool
        Constrains all taxonomic ranks in the representatives tree if set to true. 
        Only constrains the representatives to fork if false.
    localpair: bool
        If true, forces the use of the l-INS-i MAFFT algorithm instead of the default FFT-NS-2.


    Output
    ------
    chunks/unaligned/[...].fasta
        FASTA file of unaligned sequence data per chunk.
    discarded/[...].fasta
        FASTA files with discarded sequences, grouped per reason for discarding.
    supertree/representatives_tree.tre
        High-level/overarching tree containing chunk representatives.
    """

    # Importing data
    data = UniteData(fasta_path, taxon_path, length_tolerance=length_tolerance)
    # Calculating/loading in distance matrix
    distances = PairwiseDistance(matrix_file_path, indices=data.indices)

    # Creating the chunks
    chunks, disc_unidentified = data.create_chunks(min_split_depth, max_depth=max_split_depth, max_size=max_chunk_size)
    chunk_sizes = chunk_size_report(chunks)

    # Discarding distant sequences and small chunks 
    chunks, disc_distant = distances.discard_distant_sequences(chunks, strictness=outlier_strictness)
    chunks, disc_chunksize = discard_small_chunks(chunks)
    discarded_seqs = [disc_unidentified, disc_distant, disc_chunksize]
    filenames = ["unidentified", "distant", "smallchunks"]
    data.export_discarded_seqs(discarded_seqs, filenames, dir=result_dir)
    chunk_size_report(chunks)

    # Determining chunk representatives
    chunks = distances.determine_representatives(chunks, algorithm=representatives_alg)
    data.representatives_to_fasta(chunks, dir=result_dir)
    data.export_constraint_tree(chunks, dir=result_dir, full_constraint=full_constraint)

    # Generating supertree from chunk representatives
    supertree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre', mafft_path=mafft_path, raxml_path=raxml_path, localpair=localpair, dir=result_dir, constrained=True)
    non_matching = [ids[0].split("_")[0] for ids in supertree.find_non_matching_forks()]

    # Determining outgroups
    t0 = time.time()
    print("Determining outgroups for", len(chunks), "chunks")
    for i in range(len(chunks)):
        if not chunks[i].id in non_matching:
            print("|", i)
            try:
                chunks[i].outgroup = supertree.determine_outgroup(chunks[i], data) # Based on supertree
            except IndexError:
                chunks[i].outgroup = distances.determine_outgroup(chunks[i], n=1) # Fall back on distances
    t1 = time.time()
    print("Outgroups determined in", t1-t0, "seconds")

    # Saving as fasta
    data.chunks_to_fasta(chunks, exclude=non_matching, dir=result_dir)


def build_subtrees(result_dir: str, mafft_path: str, raxml_path: str, localpair: bool):
    """Generates subtrees from chunks."""

    chunk_paths = os.listdir(result_dir + "chunks/unaligned")

    # Sorting on size
    sizes = []
    for filename in chunk_paths:
        size = int(filename.split('_')[1])
        sizes.append(size)
    index_sorted = np.argsort(sizes)

    # Alignments
    i = 1
    start = time.time()
    print("Aligning", len(index_sorted), "chunks...")
    for index in index_sorted:
        t0 = time.time()
        mafft = MafftCommandline(mafft_path, localpair=localpair, thread = -1, input=result_dir + "chunks/unaligned/"+ chunk_paths[index])
        stdout, sterr = mafft()
        with open(result_dir + "chunks/aligned/" + chunk_paths[index], "w") as handle:
            handle.write(stdout)
        t1 = time.time()
        print("|", i, "    ", t1-t0, "seconds")
        i += 1
    end = time.time()
    print("Chunks aligned in", end-start, "seconds.")

    # Tree generation
    i = 1
    start = time.time()
    print("Generating", len(index_sorted), "subtrees...")
    for index in index_sorted:
        t0 = time.time()
        raxml = RaxmlCommandline(raxml_path, threads=8, sequences=result_dir + "chunks/aligned/" + chunk_paths[index], model="GTRCAT", name=chunk_paths[index][:-6], outgroup="OUTGROUP")
        raxml()
        tree = Phylo.read('RAxML_bestTree.' + chunk_paths[index][:-6], 'newick')
        Phylo.write(tree, result_dir + 'chunks/trees/' + chunk_paths[index][:-6] + '.tre', 'newick')
        files = glob.glob('*' + chunk_paths[index][:-6])
        for f in files:
            os.remove(f)
        t1 = time.time()
        print("|", i, "    ", t1-t0, "seconds")
        i += 1
    end = time.time()
    print("All trees generated in", end-start, "seconds.")


def build_supertree(result_dir: str, fasta_path: str, taxon_path: str):
    '''Inserts subtrees in representatives tree'''

    # Generating the tree
    rep_tree = RepresentativesTree(result_dir + 'supertree/representatives_tree.tre')
    backbone = Backbone(representatives_tree=rep_tree, dir=result_dir)

    # Grouping report
    data = UniteData(fasta_path, taxon_path,length_tolerance=1000)
    backbone.grouping_report(data)


def run(args: list):
    '''Processes user input and calls upon all other functions.
    For the generation of a UNITE backbone tree.'''

    sys.setrecursionlimit(2000) # For highly nested trees

    start = time.time()
    if len(sys.argv) != 9:
        print("Please provide space-separated parameter settings in the following order:")
        print("length_tolerance min_split_depth max_split_depth max_chunk_size outlier_strictness representatives_alg full_constraint localpair")
    else:
        # Setting variables
        fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
        taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
        matrix_file_path = 'data/pw_distances_6mergoogle.npy'
        mafft_path = "bin/mafft"
        raxml_path = "raxml/raxmlHPC-PTHREADS-SSE3"
        length_tolerance = float(args[1])
        min_split_depth = int(args[2])
        max_split_depth = int(args[3])
        max_chunk_size = int(args[4])
        outlier_strictness = float(args[5])
        representatives_alg = int(args[6])
        full_constraint = bool(int(args[7]))
        localpair = bool(int(args[8]))

        # Generating folder name
        result_dir = ("results/l" + str(length_tolerance) + "_s" + str(min_split_depth) + "_" + str(max_split_depth) + "_" 
        + str(max_chunk_size) + "_o" + str(outlier_strictness) + "_a" + str(representatives_alg))
        if full_constraint:
            result_dir += "_constr"
        if localpair:
            result_dir += "_localpair"
        result_dir += "/"

        # Initializing result directory
        if not os.path.isdir(result_dir[:-1]):
            os.mkdir(result_dir[:-1])
            os.makedirs(result_dir + "chunks/aligned")
            os.makedirs(result_dir + "chunks/unaligned")
            os.makedirs(result_dir + "chunks/trees")
            os.makedirs(result_dir + "supertree")
            os.makedirs(result_dir + "discarded")
        sys.stdout = open(result_dir + "log.txt", "w")
        sys.stderr = open(result_dir + "errors.txt", "w")

        # Calling functions
        divide_data(result_dir, fasta_path, taxon_path, matrix_file_path, mafft_path, raxml_path, length_tolerance, 
        min_split_depth, max_split_depth, max_chunk_size, outlier_strictness, representatives_alg, full_constraint, localpair)
        build_subtrees(result_dir, mafft_path, raxml_path, localpair)
        build_supertree(result_dir, fasta_path, taxon_path)
    
    end = time.time()
    print("Total algorithm completion time:", end-start, "seconds.")


if __name__ == '__main__':
    run(sys.argv)