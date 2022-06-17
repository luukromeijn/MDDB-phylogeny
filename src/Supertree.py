import glob, os, time
from Bio import Phylo, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Chunks import Chunk, UniteData

def remove_outgroup(tree):
    for clade in tree.root.clades:
        if clade.name == "OUTGROUP":
            tree.prune(clade)
    return tree

class Supertree:
    '''Base class for backbone/representatives tree, containing overlapping methods'''
    
    def __init__(self):
        self.tree = None # Just for syntax

    
    def grouping_report(self, data: UniteData):
        '''Shows info about how well the taxon ranks have been grouped'''

        print("Tree has", len(self.tree.get_terminals()), "leaves.")
        tax_tree = data.sh_to_tax_tree(self.tree)
        print("| " + "Rank:".ljust(10) + "Non-matching leaves:".ljust(23) + "Unique ranks:")
        for i in range(1,7):
            self.grouping_report_per_rank(tax_tree, i)

    
    def grouping_report_per_rank(self, tree, rank: int=6):
        '''Shows info about how well a specified taxon rank has been grouped.
        The tree must contain taxonomic information.'''

        all_taxons = set() # Get unique taxonomies
        for leaf in tree.get_terminals():
            taxon = self.__name_until_rank(leaf.name, rank)
            if taxon[-12:] != "unidentified" and taxon[-8:] != "OUTGROUP":
                all_taxons.add(taxon)
        
        # Get non-matching terminals
        non_matching = self.__non_matching_leaf_names(tree.root, rank)

        ranks = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
        print("|", (ranks[rank-1] + ":").ljust(10) + str(non_matching).ljust(23) + str(len(all_taxons)))


    def __non_matching_leaf_names(self, clade, rank):
        '''Recursive call of grouping value,
        Returns number of non-matching terminals'''

        count = 0
        child_1 = clade.clades[0]
        child_2 = clade.clades[1]

        # If both non-terminal, look deeper
        if len(child_1.clades) > 0 and len(child_2.clades) > 0: 
            count += self.__non_matching_leaf_names(child_1, rank)
            count += self.__non_matching_leaf_names(child_2, rank)

        # If both terminal, check match
        elif len(child_1.clades) == 0 and len(child_2.clades) == 0:
            name_1 = self.__name_until_rank(child_1.name, rank)
            name_2 = self.__name_until_rank(child_2.name, rank)
            # Skips clades with unidentified or outgroup leaves
            if name_1 != name_2 and (name_1[-12:] != 'unidentified' and name_1[-8:] != 'OUTGROUP' and
            name_2[-12:] != 'unidentified' and name_2[-8:] != 'OUTGROUP'):
                count += 1
        
        # If one terminal, compare it to the next deeper level
        else:
            if len(child_1.clades) == 0 and len(child_2.clades) > 0:
                name_x = self.__name_until_rank(child_1.name, rank)
                count += self.__non_matching_leaf_names(child_2, rank)
                children_y = child_2.clades
            elif len(child_2.clades) == 0 and len(child_1.clades) > 0:
                name_x = self.__name_until_rank(child_2.name, rank)
                count += self.__non_matching_leaf_names(child_1, rank)
                children_y = child_1.clades
            match = False
            if name_x[-12:] == 'unidentified' or name_x[-8:] == 'OUTGROUP':
                match = True
            for subchild in children_y:
                subname = self.__name_until_rank(subchild.name, rank)
                if subname == name_x:
                    match = True
            if ((children_y[0].name is None and children_y[1].name is None) or 
            (str(children_y[0].name)[-12:] == 'unidentified' and str(children_y[1].name)[-12:] == 'unidentified')):
                match = True
            if not match:
                count += 1

        return count


    def __name_until_rank(self, name, rank):
        if name is None:
            return ''
        else:
            split_on = ['_c__', '_o__', '_f__', '_g__', '_s__', '_x__']
            return name.split(split_on[rank-1])[0]


    def show(self):
        Phylo.draw(self.tree)


class RepresentativesTree(Supertree):
    '''Generates/loads supertree based on chunk representatives,
    allows operating on this supertree.'''

    def __init__(self, path: str, mafft_path: str="bin/mafft", raxml_path: str="raxml/raxmlHPC-PTHREADS-SSE3", dir: str="", localpair: bool=False, n_representatives: int=2, constrained: bool=True):
        
        start = time.time()
        try:
            self.tree = Phylo.read(path, 'newick')
        except FileNotFoundError:

            print("Calculating tree based on representatives...")

            print("| Aligning...")
            mafft = MafftCommandline(mafft_path, localpair=localpair, input=dir+"supertree/representatives_unaligned.fasta")
            stdout, sterr = mafft()
            with open(dir + "supertree/representatives_aligned.fasta", "w") as handle:
                handle.write(stdout)

            print("| Tree generation...")
            if constrained:
                # Calling raxml
                raxml = RaxmlCommandline(raxml_path, threads=8, sequences=dir+"supertree/representatives_aligned.fasta", 
                grouping_constraint = dir + "supertree/constraint.tre", model="GTRCAT", name="representatives_tree")
            else:
                raxml = RaxmlCommandline(raxml_path, threads=8, sequences=dir+"supertree/representatives_aligned.fasta", 
                model="GTRCAT", name="representatives_tree")

            raxml()
            self.tree = Phylo.read('RAxML_bestTree.representatives_tree', 'newick')
            files = glob.glob('*representatives_tree')
            for f in files:
                os.remove(f)

            Phylo.write(self.tree, dir + 'supertree/unrooted_representatives_tree.tre', 'newick')
            self.tree.root_at_midpoint()
            Phylo.write(self.tree, dir + 'supertree/representatives_tree.tre', 'newick')
        
        self.mafft_path = mafft_path
        self.raxml_path = raxml_path
        end = time.time()
        print("Representatives tree ready in", end-start, "seconds.")


    def find_non_matching_forks(self):
        '''Finds forks in supertree that do not match.'''

        distance = 0
        terminals = self.tree.get_terminals()
        matched_by_id = [[] for i in range(len(terminals))]
        for terminal in terminals:
            if terminal.name != "OUTGROUP":
                id = int(terminal.name.split('_')[0])-1
                matched_by_id[id].append(terminal)
                preterminal = self.tree.get_path(terminal)[-2]
                distance += self.tree.distance(self.tree.root, preterminal)
        distance = distance/len(terminals)

        not_matching = []
        for same_id in matched_by_id:
            if len(same_id) > 0:
                if not self.tree.is_monophyletic(same_id):
                    not_matching.append([id.name for id in same_id])

        print("Found", len(not_matching), "non-matching forks.")
        print("Average distance to preterminals in representatives tree:", distance)
        return not_matching


    def determine_outgroup(self, chunk: Chunk, unite_data: UniteData) -> Chunk:
        '''Determines outgroup for each chunk based on the representatives tree'''

        representatives = []
        for rep in chunk.representatives:
            representatives.append(chunk.id + "_" + unite_data.sequences[rep].name)
        clade = self.tree.common_ancestor(representatives) # Get clade of representatives
        parent = self.tree.get_path(clade)[-2] # Get parent

        # Choose 'sibling' node outside of ingroup with smallest distance
        min = 1000
        for sister in parent.get_terminals(): # Get sisters
            if sister.name not in representatives: # Shouldn't be representatives themselves
                dist = self.tree.distance(clade, sister) 
                if dist < min:
                    min = dist
                    if sister.name == "OUTGROUP":
                        outgroup = "OUTGROUP"
                    else:
                        outgroup = "_".join(sister.name.split("_")[1:])

        return [unite_data.name_to_index(outgroup)]


class Backbone(Supertree):
    '''Generates/loads/allows operating on backbone tree based on a chunk division'''

    def __init__(self, path: str=None, representatives_tree: RepresentativesTree=None):

        if path is not None:
            self.tree = Phylo.read(path, 'newick')
        elif representatives_tree is not None:
            
            start = time.time()
            # Finding subtree files
            subtree_paths = os.listdir("results/chunks/trees")
            subtrees = {}
            for path in subtree_paths:
                id = path.split('_')[0]
                subtrees[id] = path

            # Parsing over representatives supertree
            tree = representatives_tree.tree
            terminals = tree.get_terminals()
            matched_by_id = [[] for i in range(len(terminals))]
            for terminal in terminals:
                id = int(terminal.name.split('_')[0])-1
                matched_by_id[id].append(terminal)
            
            for same_id in matched_by_id:
                if len(same_id) > 0: # If index exists in tree
                    id = same_id[0].name.split('_')[0]
                    if tree.is_monophyletic(same_id): # Check monophyletic
                        try: # Replace representatives with subtree
                            path = subtrees[id]
                            subtree = Phylo.read("results/chunks/trees/" + path, 'newick')
                            subtree = remove_outgroup(subtree)
                            ancestor = tree.common_ancestor(same_id)
                            new_clade = subtree.root
                            ancestor.clades = new_clade.clades
                            ancestor.name = new_clade.name
                        except KeyError: # If no subtree exists, just remove the chunk id
                            print('Subtree', id, 'not found.')
                            for sh in same_id:
                                sh.name = "_".join(sh.name.split("_")[1:])
                    else: # If not monophyletic, just remove the chunk id
                        print('Found non-monophyletic representatives.')
                        for sh in same_id:
                            sh.name = "_".join(sh.name.split("_")[1:])
            
            self.tree = tree
            end = time.time()
            print("Subtrees placed in representatives tree in", end-start, "seconds.")
            Phylo.write(tree, 'results/supertree/backbone.tre', 'newick')
        else:
            raise ValueError("Please provide a path to an existing backbone or a representatives tree to generate one from.")
        