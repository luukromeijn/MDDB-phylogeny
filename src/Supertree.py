import glob, os
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
    
    def grouping_value(self):
        '''Calculates the grouping value for a given tree.
        Describes the number of monophyletic taxonomy groups,
        divided by the number of taxonomies that are present.'''

        all_taxons = set() # Get unique taxonomies
        for leaf in self.tree.get_terminals():
            all_taxons.add(leaf.name)
        
        # Call recursive function
        count = self.__recursive_grouping_value(self.tree.root)
        if type(count) == str: # If string, taxon is unique
            count = 1
        return count/len(all_taxons)


    def __recursive_grouping_value(self, clade):
        '''Recursive call of grouping value,
        Returns number of monophyletic taxonomic groups,
        Or the name of the group if this number is 1'''

        name = ""
        count = 0
        for subclade in clade.clades: # Loop through children
            if len(subclade.clades) == 0: # If terminal
                if name != subclade.name: # And name is new
                    count += 1 # Add count
                    name = subclade.name # Save name
            else: # Go deeper in recursion if not terminal
                next = self.__recursive_grouping_value(subclade)
                if type(next) == str: # If string, taxon is unique
                    if next != name:
                        count += 1 # If name is new, add and save
                        name = next
                else: # If int, add that value
                    count += next

        if count == 1:
            return name
        else:
            return count


    def show(self):
        Phylo.draw(self.tree)


class RepresentativesTree(Supertree):
    '''Generates/loads supertree based on chunk representatives,
    allows operating on this supertree.'''

    def __init__(self, path: str, mafft_path: str="bin/mafft", raxml_path: str="raxml/raxmlHPC-PTHREADS-SSE3", dir: str="", n_representatives: int=2, constrained: bool=True):
        
        # TODO improve file handling
        try:
            self.tree = Phylo.read(path, 'newick')
        except FileNotFoundError:

            print("Calculating tree based on representatives...")

            print("| Aligning...")
            mafft = MafftCommandline(mafft_path, input=dir+"supertree/representatives_unaligned.fasta")
            stdout, sterr = mafft()
            with open(dir + "supertree/representatives_aligned.fasta", "w") as handle:
                handle.write(stdout)

            print("| Tree generation...")
            if constrained:
                # Generating the constraint file (constraining all representative pairs)
                i = 0
                constraint = "("
                for seq in SeqIO.parse(dir + "supertree/representatives_aligned.fasta", 'fasta'):
                    if i == 0:
                        constraint += "("
                    constraint += seq.id + ","
                    i += 1
                    if i == n_representatives:
                        constraint = constraint[:-1] + "),"
                        i = 0
                constraint = constraint[:-1] + ");"
                file = open(dir + "supertree/constraint.tre", 'w')
                file.write(constraint)
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

            self.tree.root_at_midpoint()
            Phylo.write(self.tree, dir + 'supertree/representatives_tree.tre', 'newick')
        
        self.mafft_path = mafft_path
        self.raxml_path = raxml_path
        print("Representatives tree ready.")


    def find_non_matching_forks(self):
        '''Finds forks in supertree that do not match.'''

        terminals = self.tree.get_terminals()
        matched_by_id = [[] for i in range(len(terminals))]
        for terminal in terminals:
            id = int(terminal.name.split('_')[0])-1
            matched_by_id[id].append(terminal)

        not_matching = []
        for same_id in matched_by_id:
            if len(same_id) > 0:
                if not self.tree.is_monophyletic(same_id):
                    not_matching.append([id.name for id in same_id])

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
                    outgroup = "_".join(sister.name.split("_")[1:])
                    
        return [unite_data.name_to_index(outgroup)]


class Backbone(Supertree):
    '''Generates/loads/allows operating on backbone tree based on a chunk division'''

    def __init__(self, path: str=None, representatives_tree: RepresentativesTree=None):

        if path is not None:
            self.tree = Phylo.read(path, 'newick')
        elif representatives_tree is not None:
            
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
            Phylo.write(tree, 'results/supertree/backbone.tre', 'newick')
        else:
            raise ValueError("Please provide a path to an existing backbone or a representatives tree to generate one from.")
        