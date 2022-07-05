'''Evaluates a distance matrix using a standardized measure.'''

from DistanceData import PairwiseDistance
from Chunks import UniteData
import matplotlib.pyplot as plt

fasta_path = 'data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta'
taxon_path = 'data/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_97_10.05.2021.txt'
matrix_file_path = 'data/pw_distances_6mergoogle.npy'

data = UniteData(fasta_path, taxon_path, length_tolerance=0.2)
distances = PairwiseDistance(matrix_file_path, indices=data.indices)

scores = []
for i in range(1,7):
    chunks, _ = data.create_chunks(i, i)
    scores.append(distances.evaluate(chunks, full_output=True))

labels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
plt.boxplot(scores, showfliers=False, labels=labels)
plt.xlabel('Chunk splitting rank')
plt.ylabel('Average distance (std from mean)')
plt.title('Ingroup distances (w-metric) for different chunk divisions')
plt.savefig('dist_matrix_evaluation_w_metric.png')
