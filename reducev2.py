import os
import numpy as np
from ase.io import read, write
from scipy.cluster.hierarchy import linkage, fcluster

def calculate_rmsd(structure1, structure2):
    """Calculates the RMSD between two structures."""
    pos1 = structure1.get_positions()
    pos2 = structure2.get_positions()
    return np.sqrt(np.mean(np.sum((pos1 - pos2) ** 2, axis=1)))

def calculate_rmsd_matrix(structures):
    """Calculates the RMSD matrix between all provided structures."""
    num_structures = len(structures)
    rmsd_matrix = np.zeros((num_structures, num_structures))

    for i in range(num_structures):
        for j in range(i + 1, num_structures):
            rmsd_val = calculate_rmsd(structures[i], structures[j])
            rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd_val

    return rmsd_matrix

def select_representative_structures(structures, num_to_select):
    """Selects a specified number of representative structures based on similarity."""
    rmsd_matrix = calculate_rmsd_matrix(structures)
    rmsd_condensed = rmsd_matrix[np.triu_indices(len(structures), 1)]
    Z = linkage(rmsd_condensed, 'average')
    clusters = fcluster(Z, num_to_select, criterion='maxclust')
    
    selected_indices = []
    for cluster_id in range(1, num_to_select + 1):
        indices_in_cluster = np.where(clusters == cluster_id)[0]
        selected_indices.append(indices_in_cluster[0])  # Selects the first one from each cluster

    selected_structures = [structures[i] for i in selected_indices]
    return selected_structures

def main():
    input_folder = input("Enter the path to the folder containing the structures: ")
    output_folder = "final_representatives"

    # Create the output fold

