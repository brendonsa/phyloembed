import tempfile
import subprocess
import pandas as pd
import os
from scripts.compute_distance_matrix_phylip import pairwise_distances, write_phylip_relaxed
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
CONDA = "/home/compute2/miniconda3/bin/conda" 

def embeddings_to_score(df:pd.DataFrame, dataset:str, distance:str = 'euclidean', tree_path = None):
    """
    If tree_path is provided the trees created will be persistant and not deleted after scoring.
    """
    X = df.values
    y = df.index
    if distance == 'raw':
        distance = df.values
    else:
        distance = pairwise_distances(X,  distance)
    with tempfile.NamedTemporaryFile(suffix=".phylip", delete=False, delete_on_close=False) as tmp:
        tmp_dist_matrix = tmp.name
        write_phylip_relaxed(distance, y,tmp_dist_matrix)
    if tree_path is None:
        with tempfile.NamedTemporaryFile(suffix=".nwk", delete=False, delete_on_close=False) as tmp:
            tree_file = tmp.name
    else:
        tree_file = tree_path
    
    subprocess.run([CONDA, "run", "-n", "phyloembed",
                    "fastme", 
                    "-i",
                    tmp_dist_matrix,
                    "-o",
                    tree_file,
                    "--nni",
                    "--spr" ])
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        tmp_path = tmp.name
        tree_ref = os.path.join('./results',dataset,'tree.nwk')
        subprocess.run(
                [CONDA, "run", "-n", "phyloembed",
                    "Rscript",
                    "scripts/compare_trees.R",
                    "--t1",
                    tree_ref,
                    "--t2",
                    tree_file,
                    "--out",
                    tmp_path,
                ],
                check=True,
            )
        df = pd.read_csv(tmp_path)
        os.remove(tmp_path)
    if tree_path is None:
        os.remove(tree_file)
    os.remove(tmp_dist_matrix)

    return df