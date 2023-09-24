""" This scripts contain various modules to obtain a confidence measure
    for a given Cocktail prediction.

    There are 2 main csv inputs:
        1. "group lassoed tf df"
            This contains information of selected TF groups at various values of
            lambda (Group Lasso results) from the complete data set.
            The Cocktail prediction is present here.
            
            e.g. "BTN_group_lassoed_tf_df.csv"
                  Line of file: "lambda_2","Neurogenin,COE,HNF6"

        2. "subsampled tf df"
            This contains Group Lasso results from multiple subsamples ("runs").
        
            e.g. "BTN_combined_subsampled_tf_df.csv"
                  Line: "run_1_lambda_3","COE,HNF6"

    Gets a percentage/confidence for the overall prediction by counting the number
    of times that a TF in the predicted Cocktail appears in a subsample's inferred Cocktail.
    
    Notes:
        1. The user decides how many TFs should belong in a Cocktail ("ideal_set_size").
        We have chose 3 throughout this work.
        2. TF_set is the same as a subsample's inferred Cocktail.

    Author: Prakriti Paul Chacha
    Written: 1/20/23.

"""
import sys
import csv
from pprint import pprint

def get_lambda_and_cocktail_tfs(group_lassoed_tf_df_csv, num_tfs):
    """ Returns the first lambda that renders a set of TFs
        of size num_tfs and the names of those TFs.

        Args:
            group_lassoed_tf_df_csv: as previously described.
                Example line: "lambda_2","Neurogenin,COE,HNF6"
            num_TFs: size of potential Cocktail.

        Returns:
            lmdbda_and_cocktail_tfs: list with 2 entries:
                lambda_* and list of selected TF groups.
                e.g. ['lambda_2', ['Neurogenin', 'COE', 'HNF6']]
            None if no lambdas correspond to a set of size num_TFs.
    """
    lmdbda_and_cocktail_tfs = list()
    with open(group_lassoed_tf_df_csv, newline="") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for line in csv_reader:
            lmbda, TF_string = line[0], line[1]
            TF_string_list = TF_string.split(",")
            # In the case when no sets were rendered by lambda
            # Often at lambda_0.
            if TF_string_list == [""]:
                continue
            if len(TF_string_list) == num_tfs:
                lmdbda_and_cocktail_tfs = [lmbda, TF_string_list]
                return lmdbda_and_cocktail_tfs
        return None

def get_run(run_lambda_string):
        """ Gets the run number from run_lambda_string.
        
        Args:
            run_lambda_string.
                e.g. run_1_lambda_1
        Returns:
            run_1
        """
        return run_lambda_string.split("_lambda")[0]

def get_tf_set(tf_string):
    """ Returns a set with TFs present in tf_string.

        Args:
            tf_string.
                e.g. 'Neurogenin,COE,HNF6'
        Returns:
            {'HNF6', 'Neurogenin', 'COE'}
    """
    return set(tf_string.split(","))

def update_cocktail_confidence_dict(cocktail_confidence_dict, cocktail_tfs, TF_set):
    """ If a TF in TF_set belongs to cocktail_tfs, its count is incremented in
        cocktail_confidence_dict. Otherwise, it becomes a new key with count 1.

        Args:
            cocktail_confidence_dict: dictionary with shape {TF: counts}.
                e.g. {'Neurogenin': 46, 'HNF6': 46, 'COE': 46}
            cocktail_tfs: self-explanatory.
            TF_set: TFs in a subsample's inferred Cocktail.
        
        Returns:
            cocktail_confidence_dict: updated as described.
    """
    for tf in TF_set:
        if tf in cocktail_tfs:
            cocktail_confidence_dict[tf] += 1
        elif tf in cocktail_confidence_dict:
            cocktail_confidence_dict[tf] += 1
        else:
            cocktail_confidence_dict[tf] = 1
    return cocktail_confidence_dict

def make_lambda_and_cocktail_TF_dict(group_lassoed_tf_df_csv, TF_set_sizes):
    """ Takes a group_lassoed_tf_df_csv and makes lambda_and_cocktail_TF_dict,
        in which keys = lambda_* and values = set of selected TF Groups.

        Args:
            group_lassoed_tf_df_csv: csv file that contains information
                about selected TF groups at various lambdas.
                e.g. "BTN_group_lassoed_tf_df.csv"
                    will have a line such as "lambda_2","Neurogenin,COE,HNF6"
            
            cocktail_sizes: list of various sizes of TF Groups selected.
                e.g. [3, 4, 5]

        Returns:
            lambda_and_cocktail_TF_dict: Dictionary as described above.
                e.g. {3: ['lambda_2', ['Neurogenin', 'COE', 'HNF6']],
                      4: ['lambda_4', ['Orphan bHLH-1', 'Neurogenin', 'COE', 'HNF6']],
                      5: None...}

        Use Case:
            This dictionary is made only for "group lassoed tf df"
            (not for "subsampled tf df"). It helps the user decide what the Cocktail
            TFs are when they pick an ideal_set_size.
            i.e. cocktail_tfs = set(lambda_and_cocktail_tfs_dict[ideal_set_size][1])

            Note: The user can choose any ideal_set_size. 3 is chosen throughout this work.
    """
    # Get the lambdas required to get certain number of cocktail TFs.
    lambda_and_cocktail_tfs_dict = dict()
    for TF_set_size in TF_set_sizes:
        lmdbda_and_cocktail_tfs = get_lambda_and_cocktail_tfs(group_lassoed_tf_df_csv, TF_set_size)
        lambda_and_cocktail_tfs_dict[TF_set_size] = lmdbda_and_cocktail_tfs

    return lambda_and_cocktail_tfs_dict


def subsampling_prediction_confidence(group_lassoed_tf_df_csv, subsampled_tf_df_csv, tf_set_sizes, ideal_set_size):
    """ Performs prediction confidence routine. 
        
        Args:
            group_lassoed_tf_df_csv: as previously described.
            subsampled_tf_df_csv: as previously described.
            tf_set_sizes: list of sizes of selected TF groups (ints).
                This will be used to make "lambda_and_cocktail_tfs_dict"
                (see above for description).
            ideal_set_size: corresponds to size of Cocktail (int).
                We choose 3 throughout.

        Routine:
            Iterates through all subsamples/runs in "subsampled_tf_df_csv",
            finds sets of "ideal_set_size", keeps track of the TFs within them,
            compares these TFs with predicted Cocktail TFs
            (as found in "group_lassoed_tf_df_csv"), and updates
            "cocktail_confidence_dict" accordingly (see "update_cocktail_confidence_dict").

        Returns:
            normalized_cocktail_confidence_dict:
                dictionary with shape {TF: normalized count}
                Where normalized count is the number of times a TF in a subsample's
                inferred Cocktail is present in the predicted Cocktail divided by
                the number of subsamples/runs.
    """
    lambda_and_cocktail_tfs_dict = make_lambda_and_cocktail_TF_dict(group_lassoed_tf_df_csv, tf_set_sizes)
    
    # Keeps track of subsamples with a TF set of the ideal size.
    # Counts will be normalized with this number.
    ideal_seen_runs = set()
    seen_runs = set()

    cocktail_tfs = set(lambda_and_cocktail_tfs_dict[ideal_set_size][1])
    cocktail_confidence_dict = {cocktail_tf: 0 for cocktail_tf in cocktail_tfs}

    with open(subsampled_tf_df_csv, newline="") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for line in csv_reader:
            run, tf_set = get_run(line[0]), get_tf_set(line[1])
            seen_runs.add(run)
            if run not in ideal_seen_runs:
                if len(tf_set) == ideal_set_size:
                    ideal_seen_runs.add(run)
                    cocktail_confidence_dict = update_cocktail_confidence_dict(cocktail_confidence_dict, cocktail_tfs, tf_set)

    normalized_cocktail_confidence_dict = {TF: round(count/len(ideal_seen_runs)*100) for TF, count in cocktail_confidence_dict.items()}
    pprint(normalized_cocktail_confidence_dict)

    return normalized_cocktail_confidence_dict
