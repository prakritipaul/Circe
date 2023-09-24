""" This script contains various modules that allow for the identification
    of Pseudostates along a cell segment of the Virtual Lineage Tree given
    the Prinicipal Component and Pseudotime Information of those cells.

    Pseudostates are approximations of true Precursor, Precursor Sibling(s)
    and Mother states. They will be used to calculate Differentially Expressed
    Genes (both Precursor Identity Genes (Precursor DEGs) and DEG TFs) for the Group Lasso Model.

    Please refer to Methods of paper for Algorithm Details.
    
    Pseudocode:
        Construct an undirected weighted graph as specificed in Methods.
        Nodes are cells, edges exist between 2 cells if they are within "k"
        neighbors of each other and edge weights are Euclidean distance between
        their PC coordinates.
    
        For a k:
            For iteration 1:50:
                Do Leiden clustering
                Get cutoff index
            Take mean of cutoff indices.
        Take median of means from all k’s. This will be the final cutoff.

    Author: Prakriti Paul Chacha
    Written: 1/20/23.
"""
import csv
import decimal
import math as m
import igraph as ig
import statistics as stats
from pprint import pprint
from scipy.stats import spearmanr
from scipy.spatial import distance
from collections import OrderedDict, defaultdict

#################################################################################################################
# Helper Functions
#################################################################################################################

def make_pc_pseudotime_dict(df_file):
    """ Makes a dictionary that holds PC coordinates and pseudotimes for all cells 
        from a given df_file. which was generated from cocktail_pipeline code.

        Args:
            df_file: csv file generated from data frame in cocktail_pipeline_code
                with columns = cell_barcode, pseudotime, and PC coordinates.
            e.g. "BTN_precursor_segment_PC_pseudotime_df.csv"
        
        Returns:
            pc_pseudotime_dict:
                Ordered dictionary with keys = cell barcodes and
                values = list with pseudotime and list of PC coordinates.

            e.g. OrderedDict([('ITB.1_CAGCGACGTACGCACC',
                             [0.348199729897684,
                             [-3.3525444428361,
                              5.28799466490279...])
    """
    pc_pseudotime_dict = OrderedDict()
    with open(df_file) as df_file:
        df_reader = csv.reader(df_file, delimiter="\t")
        next(df_reader) # Drop header
        for elem in df_reader:
            cell_barcode, pseudotime, pc_coordinates = elem[0], float(elem[1]), list(map(float, elem[2:]))
            pc_pseudotime_dict[cell_barcode] = [pseudotime, pc_coordinates]
    return pc_pseudotime_dict

def get_all_pc_vectors(pc_pseudotime_dict):
        """ Extracts pc_vectors from all cells in a pc_pseudotime_dict.
            Returns a list of pc_vectors (list of floats).
        """
        dict_items = list(pc_pseudotime_dict.items())
        pc_vectors = [dict_item[1][1] for dict_item in dict_items]
        return pc_vectors

def get_distance_matrix(pc_vec_list, k):
    """ Given a list of PC vectors and k, function gets Euclidean
        distance between a given cell and up to k of its neighbors
        to generate a distance_matrix.
        
        Args:
            pc_vec_list:
                List of PC vectors (Each PC vector is a list of floats).
            k: number of neighbors (int).
                We will calculate the Euclidean distance between a cell and upto k
                of its neighboring cells.

        Returns:
            distance_matrix: list of lists.
                Each entry is a pair-wise Euclidean distance between cells.
    """
    # Initialize n x n square matrix where n = number of cells.
    distance_matrix = [[0.0] * len(pc_vec_list) for _ in range(len(pc_vec_list))]
    for i in range(len(pc_vec_list)):
        for j in iter(num for num in range(max(0, i-k), min(i+1+k, len(pc_vec_list))) if num != i):
            if distance_matrix[i][j] == 0.0:
                pairwise_distance = distance.euclidean(pc_vec_list[i], pc_vec_list[j])
                distance_matrix[i][j] = pairwise_distance
            else:
                distance_matrix[j][i] = distance_matrix[i][j]
    return distance_matrix

def make_all_edges_list(distance_matrix):
    """ Iterates through distance_matrix and returns edge_list, a list that
        contains (i, j, weight) tuples.
    """
    edge_list = []
    for i in range(len(distance_matrix)):
        for j in range(len(distance_matrix)):
            edge_list_entry = (i, j, distance_matrix[i][j])
            edge_list.append(edge_list_entry)
    return edge_list

def make_weighted_undirected_G_from_tuple_list(edge_list):
    """ Makes a weighted undirected graph given edge_list.
    """
    G = ig.Graph()
    G = ig.Graph().TupleList(edge_list, directed=False, weights=True)
    return G

def do_leiden_clustering(weighted_undirected_G, resolution_parameter=1, n_iterations=10):
    """ Performs Leiden Clustering on a weighted undirected graph by optimizing
        the modularity measure. 

        Args:
            weighted_undirected_G: an iGraph Graph object.
            resolution_parameter: self-explanatory.
            n_iterations: the number of iterations to iterate the Leiden algorithm.

        Returns: 
            clustering: a VertexClustering object.
    """
    clustering = weighted_undirected_G.community_leiden(objective_function='modularity', weights='weight', resolution_parameter=resolution_parameter, n_iterations=n_iterations)
    return(clustering)

def get_cutoff_index(cluster_label_list, correct_label):
    """ Given cluster_label_list (a list of 0's and 1's) and correct_label (0),
        finds the index at which the label changes from correct_label (0)
        to another label (1).

        Use Case: The 0's label Pseudostate cells and 1's label non-Pseudostate cells.
            This index demarcates the end of the Pseudostate along the cell segment.
    """
    for (index, label) in enumerate(cluster_label_list):
        if label != correct_label:
            return index

def get_num_disagreements(cluster_label_list, verbose):
    """ Counts the number of times cluster label for cell i differs from cell i+1.
        Function returns this value.
        
        Args:
            cluster_label_list: list of 0's and 1's.
            verbose: boolean. Determines whether to print number of clusters and 
                number of disagreements or not.
        Returns:
            nunber of diagreements: As described.

        Use Case:
            The best clustering will have the lowest number of disagreements. This implies
            the longest streak of the same cluster labels and thus greatest homogeneity
            within the clusters. 

            The number of disagreements is used as a tie-breaker between two clusterings
            with the same modularity.
    """
    assert len(set(cluster_label_list)) > 0, "Expected more than 1 cluster"

    num_disagreements = 0
    for i in range(len(cluster_label_list)-1):
        if cluster_label_list[i] != cluster_label_list[i+1]:
            num_disagreements += 1

    if verbose: print(f"clusters = {clusters} and num_disagreements = {num_disagreements}")
    return(num_disagreements)

def get_top_cutoff_tuple(G, n_leiden_iterations=10, verbose=False):
    """ Gets the best/top cutoff index of a pseudotimed population based on a
        clustering's modularity and disagreement ratio (DR).

        Args:
            G: a weighted undirected graph.
            n_leiden_iterations: number of iterations to be used by Leiden algorithm.
            verbose: boolean. Determines whether to print various outputs
                of function or not.

        Routine:
            1. First performs Leiden clustering with lowest resolution.
            2. Increments resolution by steps of 0.01 until there are 3 clusters.
            3. While there are 2 clusters, gets the...
                a. Cutoff index (marks the start of the next cluster,
                   corresponds to the non-Pseudostate).
                b. Modularity.
                c. Number of disagreements
                
                Saves (resolution, cutoff, number of disagreements, modularity)
                in output_tuple_list.

        Returns:
            top_cutoff_tuple: tuple that contains (resolution_parameter, cutoff_index,
                num_disagreements, modularity) of the top cutoff_index. The top
                cutoff_index corresponds to the best clustering based on the highest
                modularity and fewest number of disagreements.
    """
    output_tuple_list = list()

    resolution_parameter = 0.01
    num_clusters = len(do_leiden_clustering(G, resolution_parameter, n_leiden_iterations))
    if verbose: print(f"Initial num_clusters = {num_clusters} with resolution_parameter = {resolution_parameter}")

    while num_clusters != 3:
        if verbose: print(f"\nStart of while loop: num_clusters = {num_clusters}")
        resolution_parameter += 0.01
        clustering = do_leiden_clustering(G, resolution_parameter=resolution_parameter, n_iterations=n_leiden_iterations)
        num_clusters, modularity, membership = len(clustering), clustering.modularity, clustering.membership

        if verbose: print(f"next resolution_parameter = {resolution_parameter}")
        if verbose: print(f"num_clusters post-clustering = {num_clusters}")
        if num_clusters == 2:
            # Get components of output tuple: Correct Label, Cutoff, and DR.
            correct_label = 0
            cutoff_index = get_cutoff_index(membership, correct_label)
            num_disagreements = get_num_disagreements(membership, False)
            
            output_tuple = (resolution_parameter, cutoff_index, num_disagreements, modularity)
            if verbose: print(f"Entered if: resolution_parameter = {resolution_parameter}, cutoff_index = {cutoff_index}, num_disagreements = {num_disagreements}, modularity = {modularity}")
            output_tuple_list.append(output_tuple)
    # Get cutoff that corresponds to the highest modularity and lowest num_disagreements
    top_cutoff_tuple = max(output_tuple_list, key=lambda x: (x[3], -x[2]))
    return top_cutoff_tuple

#################################################################################################################
# Main Pipeline Function
#################################################################################################################

def get_pseudostates_pipeline(df_file, n_leiden_iterations, n_routine_iterations, verbose=False):
    """ Given pseudotime and PC coordinates of cells along a segment, finds pseudostates
        using Leiden clustering.

        Args:
            df_file: csv file generated from data frame in cocktail_pipeline_code
                with columns = cell_barcode, pseudotime, and PC coordinates.
                e.g. "BTN_precursor_segment_PC_pseudotime_df.csv"
            n_leiden_iterations: number of iterations to be used by Leiden algorithm.
            n_routine_iterations: number of times to iterate the Leiden algorithm
                for a given k.
                Note: k = 5, 10, and 15% of pseudotimed cells.
            verbose: boolean. Determines whether to print various outputs
                of function or not.
        
        Routine:
            Refer to docstring of "get_top_cutoff_tuple".

            For a k, performs "get_top_cutoff_tuple"  n_routine_iterations times and
            then takes the median of the resultant cutoff indices. It does this
            for the 3 aforementioned k's (final_cutoffs) and gets the average (final_mean_cutoff).

        Returns:
            final_cutoffs: list of 3 ints.
            final_mean_cutoff: int.
        
    """
    pseudotime_dict = make_pc_pseudotime_dict(df_file)
    all_pc_vectors = get_all_pc_vectors(pseudotime_dict)
    distance_function = distance.euclidean
    k_list = [round(0.05*len(pseudotime_dict)), round(0.1*len(pseudotime_dict)), round(0.15*len(pseudotime_dict))]
    
    final_cutoff_tuples = list()
    final_cutoffs = list()

    for k in k_list:
        if verbose: print(f"###### Doing for k = {k}\n\n")
        distance_matrix = get_distance_matrix(all_pc_vectors, k)
        edge_list = make_all_edges_list(distance_matrix)
        g = make_weighted_undirected_G_from_tuple_list(edge_list)

        # Contains all the info.
        top_routine_cutoff_tuples = list()
        # Contains just the number.
        top_routine_cutoffs = list()

        for i in range(n_routine_iterations):
            # REMOVE
            # print(f"Doing for k = {k}, i = {i}")
            # Get a top cutoff for each iteration.
            if verbose: print(f"######## Doing for i = {i} ########")
            top_routine_cutoff_tuple = get_top_cutoff_tuple(g, n_leiden_iterations=n_leiden_iterations, verbose=verbose)
            top_routine_cutoff = top_routine_cutoff_tuple[1]

            top_routine_cutoff_tuples.append(top_routine_cutoff_tuple)
            top_routine_cutoffs.append(top_routine_cutoff)

        final_cutoff = stats.median(top_routine_cutoffs)
        if verbose: print(f"top_routine_cutoff = {top_routine_cutoff}")
        final_cutoff_tuples.append(top_routine_cutoff_tuples)
        final_cutoffs.append(final_cutoff)

    # Take the mean across the 3 k's.
    final_mean_cutoff = round(stats.mean(final_cutoffs))

    if verbose: print(f"\nfinal cutoff answer = {final_mean_cutoff}")
    if verbose: print("\nfinal_cutoffs")
    if verbose: print(final_cutoffs)
    if verbose: print("\nfinal_cutoff_tuples")
    if verbose: pprint(final_cutoff_tuples)

    return final_cutoffs, final_mean_cutoff
