import glob, os
from Bio import Phylo
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Chunks import Chunk, UniteData

def remove_outgroup(tree):
    for clade in tree.root.clades:
        if clade.name == "OUTGROUP":
            tree.prune(clade)
    return tree


class RepresentativesTree:
    '''Generates/loads supertree based on chunk representatives,
    allows operating on this supertree.'''

    def __init__(self, path: str, mafft_path: str="bin/mafft", raxml_path: str="raxml/raxmlHPC-PTHREADS-SSE3"):
        
        # TODO improve the path handling and everything related to that
        try:
            self.tree = Phylo.read(path, 'newick')
        except FileNotFoundError:

            print("Calculating tree based on representatives...")

            print("| Aligning...")
            mafft = MafftCommandline(mafft_path, input="results/supertree/representatives_unaligned.fasta")
            stdout, sterr = mafft()
            with open("results/supertree/representatives_aligned.fasta", "w") as handle:
                handle.write(stdout)

            print("| Tree generation...")
            raxml = RaxmlCommandline(raxml_path, threads=8, sequences="results/supertree/representatives_aligned.fasta", model="GTRCAT", name="representatives_tree")
            raxml()
            self.tree = Phylo.read('RAxML_bestTree.representatives_tree', 'newick')
            Phylo.write(self.tree, 'results/supertree/representatives_tree.tre', 'newick')
            files = glob.glob('*representatives_tree')
            for f in files:
                os.remove(f)
        
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
                    # for id in same_id:
                    #     id.color = 'red'
                    not_matching.append([id.name for id in same_id])

        return not_matching


    def fix_non_matching_forks(self): # TODO not working atm
        '''Regenerates the representatives tree, 
        placing a grouping constraint on the non-matching forks.'''

        not_matching = self.find_non_matching_forks()

        if len(not_matching) > 0:
            constraints = ""
            for same_id in not_matching:
                new_constraint = "("
                for sh in same_id:
                    new_constraint += sh + ","
                constraints += "," + new_constraint[:-1] + ")"
            constraints = "(" + constraints[1:] + ");"
            file = open("results/supertree/constraint.tre", "w")
            file.write(constraints)

            print("Regenerating supertree, constraining", len(not_matching), "non-matching forks...")

            raxml = RaxmlCommandline(self.raxml_path, threads=8, 
            sequences="results/supertree/representatives_aligned.fasta", 
            model="GTRCAT", grouping_constraint="results/supertree/constraint.tre", 
            name="constrained_representatives_tree")

            raxml() #TODO There is an error here
            
            tree = Phylo.read('RAxML_bestTree.' + "constrained_representatives_tree", 'newick')
            Phylo.write(tree, 'results/supertree/' + "constrained_representatives_tree.tre", 'newick')
            files = glob.glob('*' + "constrained_representatives_tree")
            for f in files:
                os.remove(f)

            not_matching = self.find_non_matching_forks()
            if len(not_matching) > 0:
                raise RuntimeError("Couldn't fix", len(not_matching), "forks.")
            else:
                print("Regeneration with constraints done.")


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
                    
        return [unite_data.get_seq_index(outgroup)]


class Backbone:
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
                        print('Representatives', id, 'not monophyletic.')
                        for sh in same_id:
                            sh.name = "_".join(sh.name.split("_")[1:])
            
            self.tree = tree
            Phylo.write(tree, 'results/supertree/backbone.tre', 'newick')
        else:
            raise ValueError("Please provide a path to an existing backbone or a representatives tree to generate one from.")

    
    def show(self):
        Phylo.draw(self.tree)
        