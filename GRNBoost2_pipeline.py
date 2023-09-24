""" This script contains various modules that allow us to infer Cocktails
    using a nonlinear Gene Regulatory Network inference algorithm
    GRNBoost2 [1].

    The purpose of these modules is to compare Cocktail predictions derived
    from nonlinear models with those from the linear Group Lasso model, as TF-gene
    relationships may be better modeled as nonlinear ones.

    Please refer to Methods of paper to better understand how the GRNBoost2 algorithm works.
    The key thing to remember is that the algorithm outputs an importance measure
    for every TF-gene pair, which is the weight that the TF has in the prediction
    of the corresponding gene.

    The following analyses are accomplished by the various modules:
        1. For each Candidate Cocktail TF, its importance scores with all other
           genes are summed, rendering a “Total Importance” Score and the TFs
           are subsequently ranked based on these values.
           The top 3 Candidate Cocktail TFs are considered the Cocktail.

        2. For each DEG, the top 3 TFs with the highest importance scores are noted.
           The number of times that a Candidate Cocktail TF appeared in a top 3
           list of a DEG is counted, normalized with the number of DEGs, and is used
           to rank TFs by this percentage. This is the "Percent Covered" Score.
           Once again, the top 3 Candidate Cocktail TFs are considered the Cocktail.

    Finally, as GRNBoost is not a deterministic algorithm, it was run 50 times
    on the same dataset per cell type and a standard deviation for the values
    above is obtained.

    The ensuing results demonstrate good agreement between the two Cocktail
    predictions, which builds computational evidence for the correctness of the
    Group Lasso model and the strengh of the assumptions used in this work.

    [1] Moerman, Thomas, et al. "GRNBoost2 and Arboreto: efficient and scalable
        inference of gene regulatory networks." Bioinformatics 35.12 (2019): 2159-2161
    
    Author: Prakriti Paul Chacha
    Written: 1/20/23.
"""
import os
import sys
import csv
import glob
import pickle
import math as m
import pandas as pd
import numpy as np
from pprint import pprint
from collections import OrderedDict, defaultdict

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

#################################################################################################################
# Helper Functions
#################################################################################################################

# TF_tsv = "./helper_files/Ciona_khid_TF.tsv"
def make_munging_dict(khid_tf_tsv_file):
    """ Takes in tsv file containing khids and TF Ortholog Names and creates
        munging dict with shape {khid: TF ortholog}.

        Munging dict will be used to convert khids to Ortholog Names in other
        dictionaries such as total_tf_importance_dict and sorted_DEG_regulator_dict
        for each of interpretation.
    """
    munging_dict = dict()
    with open(khid_tf_tsv_file) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        for line in tsv_reader:
            khid, TF_ortholog = line[0], line[1]
            munging_dict[khid] = TF_ortholog
    return munging_dict

def do_GRNBoost2(expression_count_matrix_tsv, tf_khid_csv, outname, min_i, max_i):
    """ Function runs GRNBoost2 n (=max_i-min_i+1) times and creates n GRNs,
        which are then exported as tsv files with name outname_*.

        Args:
            expression_count_matrix_tsv: tsv file generated from R data frame
                with rows = cell barcodes and columns = DEGs + Candidate Cocktail TFs.
                Entries are counts from single-cell counts matrix.
            tf_khid_csv: csv file that contains khids of Candidate Cocktail TFs.
            outname: out directory of GRN tsv names.
                e.g. ./GRN/BTN/BTN -> "./GRN/BTN/BTN_1.tsv"
            min_i = First value of i
            max_i = Last value of i

        Returns:
            n number of GRN tsv files.

        Use Case:
            These GRNs are processed in downstream functions to accomplish
            aforementioned analyses to infer Cocktails.
    """
    for i in range(min_i, max_i+1):
        print(f"Building GRN #{i}")
        # Get expression matrices.
        expression_mat = pd.read_csv(expression_count_matrix_tsv, sep="\t")

        # Get TF names.
        tf_khid_names = load_tf_names(tf_khid_csv)

        # Build the GRN.
        grn = grnboost2(expression_data=expression_mat,
                        tf_names=tf_khid_names)
        # Write it out.
        grn_outname = outname + f"_{i}.tsv"
        grn.to_csv(grn_outname, sep='\t', index=False, header=False)

def make_grn_dict(grn_file):
    """ Given a grn_file (output of GRNBoost2) makes a dictionary with shape:
        {(TF, target): [Importance Score]}

        Args:
            grn_file: tsv file that contains GRN output of GRNBoost2.
                Note: Each line of tsv file consists of TF (Candidate Cocktail TF),
                    TF's target (Precursor DEG), and Importance Score.

        Returns:
            grn_dict: As desribed. Key = tuple with TF and target;
                Value = List containing Importance Score.
        
        Use Case:
            grn_dicts are created from each tsv output of do_GRNBoost2, which are
            then subsequently merged into 1 dictionary (master_dict), which
            will thus have information from all runs of GRNBoost2.

    """
    grn_dict = defaultdict(list)

    with open(grn_file, newline="") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter="\t")
        for row in csv_reader:
            tf_target, importance = (row[0], row[1]), float(row[2])
            grn_dict[tf_target].append(importance)

    return grn_dict

def merge_dicts(master_dict, grn_dict):
    """ Merges a master_dict with grn_dict to create an updated master_dict. 
        master_dict is a dict that is already merged with at least 2 other grn_dicts.

        Args:
            master_dict, grn_dict:
                Dictionaries of shape: {(TF, target): [Importance Score(s)]}
                Note: While grn_dict will only have 1 Importance Score,
                    master_dict will have 1 or more Importance Scores.
        
        Routine: Merges the dictionaries.
        
        Returns:
            master_dict: merged dictionary/updated version of input master_dict.

        Note: If both dictionaries have the same key (tf, target), then function
            appends value (importance) of grn_dict to value of MASTER_dict.
            Otherwise, function makes a new key:value pair.

        e.g.
            Inputs:
                master_dict = {('KH2012:KH.C14.377', 'KH2012:KH.C12.204'): [87.88612638],
                               ('KH2012:KH.C6.129', 'KH2012:KH.C12.204'): [71.02263476]}

                grn_dict =    {('KH2012:KH.C14.377', 'KH2012:KH.C12.204'): [90.49058166],
                               ('KH2012:KH.C6.129', 'KH2012:KH.C12.204'): [56.03677012],
                               ('KH2012:KH.C2.42', 'KH2012:KH.S618.6'): [32.84915809]}
            Output:
                master_dict = {('KH2012:KH.C14.377', 'KH2012:KH.C12.204'): [87.88612638, 90.49058166]
                               ('KH2012:KH.C6.129', 'KH2012:KH.C12.204'): [71.02263476, 56.03677012],
                               'KH2012:KH.C2.42', 'KH2012:KH.S618.6'): [32.84915809]}

    """
    for tf_target in grn_dict.keys():
        if tf_target in master_dict:
            master_dict[tf_target] = master_dict[tf_target] + grn_dict[tf_target]
        else:
            master_dict[tf_target] = grn_dict[tf_target]
    return master_dict

def make_master_dict_stats(master_dict):
    """ Makes master_dict_stats, which contains the mean and standard devation
        of Importance Scores for every TF-target pair.

        Args:
            master_dict: Dictionary with shape {(TF, target): [Importance Score(s)]}

        Returns:
            master_dict_stats: As described.
                It has shape {(TF, DEG): (u(I), sd(I))}, where u = mean,
                sd = standard deviation, and I is Importance Score.
                Both key and values are tuples.
    """
    master_dict_stats = defaultdict()
    
    for tf_target in master_dict:
        importances = master_dict[tf_target]
        mean, stdev = np.mean(importances), np.std(importances)
        master_dict_stats[tf_target] = (mean, stdev)
    return master_dict_stats

def make_DEG_regulator_dict(master_dict_stats):
    """ Given master_dict_stats, makes DEG_regulator_dict, which contains
        information about how each Precursor DEG is regulated by each Candidate
        Cocktail TF (mean and standard deviation of Importance Scores).

        Args:
            master_dict_stats: as previously described.
                It has shape {(TF, DEG): (u(I), sd(I))}

        Returns: DEG_regulator_dict:
            It has shape {DEG: [(TF_1, u(I), sd(I)), (TF_2, u(I), sd(I))...]}
            for all DEGs and all TFs. 
            key = TF and value = list of Importance Score tuples.

        Use Case:
            The importance score tuples will be sorted based on the mean
            to decide which TFs are the most important in regulating the DEG.
    """
    DEG_regulator_dict = defaultdict(list)
    for (TF, DEG) in master_dict_stats:
        importance_mean, importance_sd = master_dict_stats[(TF, DEG)][0], master_dict_stats[(TF, DEG)][1]
        DEG_regulator_dict[DEG].append((TF, importance_mean, importance_sd))
    return DEG_regulator_dict

def make_sorted_DEG_regulator_dict(DEG_regulator_dict):
    """ Given a DEG_regulator_dict, sorts the TF regulators for each DEG based
        on their Importance Score means.

        Args:
            DEG_regulator_dict: as previously described.
                It has shape {DEG: [(TF_1, u(I), sd(I)), (TF_2, u(I), sd(I))...]}

        Returns:
            sorted_DEG_regulator_dict:
                It has the same shape as above.
                The values have the Importance Score tuples sorted in nonincreasing
                order based on the means.

        Use Case:
            This dictionary will help us find the top n TF regulators for each DEG.
    """
    sorted_DEG_regulator_dict = dict()
    for DEG in DEG_regulator_dict:
        importance_tuples = DEG_regulator_dict[DEG]
        sorted_importance_tuples = sorted(importance_tuples, key=lambda importance_tup:importance_tup[1], reverse=True)
        sorted_DEG_regulator_dict[DEG] = sorted_importance_tuples
    return sorted_DEG_regulator_dict

def make_top_n_regulator_dict(sorted_DEG_regulator_dict, n):
    """ Takes the first n TF regulators for each DEG from sorted_DEG_regulator_dict.
        
        Args:
            sorted_DEG_regulator_dict: as described above.
            It has shape {DEG: [(TF_1, u(I), sd(I)), (TF_2, u(I), sd(I))...]}

            n: number of top TF regulators wanted.
        
        Returns:
            top_n_regulator_dict:
            It has shape {DEG: [TF_1, TF_2... TF_n]}
            key = DEG and value = list of n top regulators.
    """
    top_n_regulator_dict = defaultdict(list)
    for DEG in sorted_DEG_regulator_dict:
        top_n_regulators = [TF for (TF, importance_mean, importance_sd) in sorted_DEG_regulator_dict[DEG][:n]]
        top_n_regulator_dict[DEG] = top_n_regulators
    return top_n_regulator_dict

def make_total_TF_importance_dict(master_dict_stats):
    """ For each TF, gets its "Total Importance" Score, i.e. sum of its Importance Means
        from all all DEGs and the combined standard deviation.

        Args:
            master_dict_stats: as previously described.
                It has shape {(TF, DEG): (u(I), sd(I))}

        Returns:
            total_TF_importance_dict:
                It has shape: {TF: [total_importance_mean, total_importance_sd]}
                key = TF and value = list with Total importance mean and
                standard deviation.

        Use Case:
            Help determine which the most important TFs are in the GRN.
    """
    total_TF_importance_dict = dict()
    for (TF, DEG) in master_dict_stats:
        importance_info = master_dict_stats[(TF, DEG)]
        importance_mean, importance_var = importance_info[0], (importance_info[1])**2
        if TF not in total_TF_importance_dict:
            total_TF_importance_dict[TF] = list()
            total_TF_importance_dict[TF].extend([importance_mean, importance_var])
        else:
            total_TF_importance_dict[TF][0] += importance_mean
            total_TF_importance_dict[TF][1]+= importance_var

    # Get the standard deviations.
    total_TF_importance_dict = {TF:[mean, m.sqrt(var)] for (TF, (mean, var)) in total_TF_importance_dict.items()}
    # Just show 2 significant figures.
    total_TF_importance_dict = {TF:[float("%.2f" % mean), float("%.2f" % sd)] for (TF, (mean, sd)) in total_TF_importance_dict.items()}
    return total_TF_importance_dict

def munge_sorted_DEG_regulator_dict(sorted_DEG_regulator_dict, munging_dict):
    """ Replaces TF khids with Orthologs, for ease of interpretation.

        Args:
            sorted_DEG_regulator_dict: as previously described.
                It has shape {DEG: [(TF_1, u(I), sd(I)), (TF_2, u(I), sd(I))...]}

            munging_dict:
                Dictionary with keys = khids and values = Ortholog.
                e.g. BTN_TF_KHID_TO_GHOST_DICT = {"KH2012:KH.C7.269": "Orphan bHLH-1",
                                                  "KH2012:KH.C2.957": "msxb"...}
        Routine:
            Iterates through all TF khids in sorted_DEG_regulator_dict and
            replaces them with their Orthologs using munging_dict.

        Returns:
            munged_sorted_DEG_regulator_dict:
                It is the same as sorted_DEG_regulator_dict, just with Orthologs.
    """
    munged_sorted_DEG_regulator_dict = dict()
    for DEG in sorted_DEG_regulator_dict:
        munged_regulating_TF_list = list()
        regulating_TF_list = sorted_DEG_regulator_dict[DEG]
        for tup in regulating_TF_list:
            munged_tup = (munging_dict[tup[0]], tup[1], tup[2])
            munged_regulating_TF_list.append(munged_tup)
        munged_sorted_DEG_regulator_dict[DEG] = munged_regulating_TF_list
    return munged_sorted_DEG_regulator_dict

def get_percent_top(TF, top_n_regulator_dict):
    """ Returns the percent of times a given TF shows up as a
        top TF regulator for DEGs in "top_n_regulator_dict".

        Args:
            TF: self-explanatory.
            top_n_regulator_dict: as previously described.
            It has shape {DEG: [TF_1, TF_2... TF_n]}

        Returns:
            Percent of times a given TF shows up as a
            top TF regulator for DEGs in "top_n_regulator_dict", or
            "Percent Covered" Score.

    """
    percent_top = 0
    for DEG in top_n_regulator_dict:
        if TF in top_n_regulator_dict[DEG]:
            percent_top += 1
    percent_top /= len(top_n_regulator_dict)
    percent_covered = float("%.2f" % percent_top)
    return percent_covered

def get_sorted_percent_top_list(top_n_regulator_dict, munged_total_tf_importance_dict):
    """ Given a top_n_regulator_dict, sorts TFs based on Percent Covered Scores.

        Args:
            top_n_regulator_dict: as previously described.
                It has shape {DEG: [TF_1, TF_2... TF_n]}

            munged_total_tf_importance_dict:
                It is the same as total_tf_importance_dict, just with
                keys = TF Orthologs instead of TF khids.
                It has shape {TF: [total_importance_mean, total_importance_sd]}.

        Routine:
            Iterates through all DEGs in munged_total_tf_importance_dict,
            keeps count of how many times a TF appears in top_n_regulator_dict,
            then divides by total number of TFs to get a "Percent Covered" Score.

        Returns:
            sorted_percent_top_list:
                Sorted List of tuples- each tuple consists of (TF, Percent Covered Score)
                TFs are sorted by Percent Covered Score in nonincreasing order.
    """
    percent_top_list = list()
    for TF in munged_total_tf_importance_dict:
        percent_top_list.append((TF, get_percent_top(TF, top_n_regulator_dict)))
    sorted_percent_top_list = sorted(percent_top_list, key=lambda tup:tup[1], reverse=True)
    return sorted_percent_top_list

#################################################################################################################
# Main Pipeline Function
#################################################################################################################

def GRNBoost2_pipeline(grn_path, khid_tf_tsv_file, n):
    """ Performs pipeline and returns sorted_total_importance_tfs (TFs ordered by Total
        Importance Scores) and percentage covered by each TF as a Top Regulator.
    """
    grn_files = glob.glob(grn_path)

    master_dict = defaultdict(list)
    i = 1
    for grn_file in grn_files:
        munged_name = grn_file.split("_")[-1]
        grn_dict = make_grn_dict(grn_file)
        master_dict = merge_dicts(master_dict, grn_dict)
        master_dict_stats = make_master_dict_stats(master_dict)
        i += 1

    munging_dict = make_munging_dict(khid_tf_tsv_file)

    # Get Total Importance of each TF in the GRN.
    total_tf_importance_dict = make_total_TF_importance_dict(master_dict_stats)
    # Replace khids with TF Ortholog Names.
    munged_total_tf_importance_dict = {munging_dict[key]:value for (key, value) in total_tf_importance_dict.items()}
    sorted_total_importance_tfs = sorted(munged_total_tf_importance_dict.items(), key=lambda x:x[1][0], reverse=True)

    # Explore what regulates the different DEGs.
    DEG_regulator_dict = make_DEG_regulator_dict(master_dict_stats)
    sorted_DEG_regulator_dict = make_sorted_DEG_regulator_dict(DEG_regulator_dict)
    munged_sorted_DEG_regulator_dict = munge_sorted_DEG_regulator_dict(sorted_DEG_regulator_dict, munging_dict)

    # Get the Top n Regulators for each DEGs.
    top_n_regulator_dict = make_top_n_regulator_dict(munged_sorted_DEG_regulator_dict, n)
    # How many times does each TF show up as a Top Regulator?
    sorted_percent_top_list = get_sorted_percent_top_list(top_n_regulator_dict, munged_total_tf_importance_dict)

    print("sorted_total_importance_tfs")
    pprint(sorted_total_importance_tfs)

    print("\nsorted_percent_top_list")
    pprint(sorted_percent_top_list)
    return sorted_total_importance_tfs, sorted_percent_top_list
