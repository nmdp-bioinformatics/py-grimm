import os
import sys
import time
import json
from os import PathLike
from typing import Iterable, Union, Dict

import pandas as pd
from grim import grim

from grma.match import Graph as MatchingGraph
from grma.match.donors_matching import DonorsMatching, _init_results_df
from grma.match.graph_wrapper import Graph
from grma.utilities.utils import print_time, donor_mismatch_format

GRIM_DEFAULT_OUTPUT_PATH = "./output/don.pmug"
GRIM_RESULT_DIR_FIELD = "imuptation_out_path"
GRIM_RESULT_GENO_FILE_FIELD = "imputation_out_umug_freq_filename"


def run_grim(config_file_path=""):
    """"
    This function applies grim imputation with the default/user's configuration file.

    :param config_file_path: A path to configuration file for grim imputation. default is grim default config.
    :return: The path to the imputation output file.
    """
    if not config_file_path:
        grim.graph_freqs()
        grim.impute()
        return GRIM_DEFAULT_OUTPUT_PATH

    if not (os.path.exists(config_file_path) and os.path.isfile(config_file_path)):
        raise FileExistsError(f"Can't find configuration file in {config_file_path}.")

    with open(config_file_path, "r") as f:
        conf = json.load(f)

    # Finds the path of the output of grim for
    dir_path = conf[GRIM_RESULT_DIR_FIELD]
    file_name = conf[GRIM_RESULT_GENO_FILE_FIELD]
    path = os.path.join(dir_path, file_name)

    grim.graph_freqs(config_file_path)
    grim.impute(config_file_path)

    return path


def search_in_levels(patient_id: int, g_m: DonorsMatching, donors_info: Iterable, threshold: float,
                     cutof: int, classes: Iterable, subclasses: Iterable):
    """"
    This function gets the information about the patient and other settings, and search for donors in levels.
    First, it will search only for 10 matches donors (compare only genotypes). If there aren't enough donors,
    the function will search for 9 matches donors (compare classes), and then 8-7 matches (subclasses)

    :param patient_id: The id of the patient we want to search his/her matches.
    :param g_m: The patients graph with the info of the current patient (and maybe other patients).
    :param donors_info: An iterable of fields from the database to include in the results. default is None.
    :param threshold: Minimal score value for a valid match. default is 0.1.
    :param cutof: Maximum number of matches to return. default is 50.
    :param classes: An iterable with the all the possible classes of the patient.
    :param subclasses: An iterable with the all the possible subclasses of the patient.
    :return: A pandas.DataFrame with the matches for this patient.
    """

    matched = set()  # set of donors ID that have already matched for this patient
    results_df = _init_results_df(donors_info)  # initialize the df according to the given fields

    # We can give to this function the genotypes instead
    g_m.find_geno_candidates_by_genotypes(patient_id)
    matched, count, results_df = g_m.score_matches(0, results_df, donors_info, patient_id, threshold, cutof, matched)

    if len(matched) >= cutof:
        return results_df

    g_m.find_geno_candidates_by_classes(classes)
    matched, count, results_df = g_m.score_matches(1, results_df, donors_info, patient_id, threshold, cutof, matched)

    if len(matched) >= cutof:
        return results_df

    g_m.find_geno_candidates_by_subclasses(subclasses)

    # loop over possible mismatches: 2, 3.
    for mismatches in range(2, 4):
        matched, count, results_df = g_m.score_matches(mismatches, results_df, donors_info,
                                                       patient_id, threshold, cutof, matched)

    return results_df


def find_matches(imputation_filename: Union[str, PathLike], match_graph: Graph,
                 search_id: int = 1, donors_info: Iterable[str] = [],
                 threshold: float = 0.1, cutof: int = 50,
                 verbose: bool = False, save_to_csv: bool = False,
                 calculate_time: bool = False):
    """
    The main function responsible for performing the matching.
    Note: for each patient, if a donor has been found as a
    match in an early stage (0, 1, or 2 mm), he will not be searched as a match for the further mismatches.

    :param imputation_filename: Path to the output file of the imputation made by grim.
    :param match_graph: A Graph object from grma.match
    :param search_id: An integer identification of the search. default is 0.
    :param donors_info: An iterable of fields from the database to include in the results. default is None.
    :param threshold: Minimal score value for a valid match. default is 0.1.
    :param cutof: Maximum number of matches to return. default is 50.
    :param verbose: A boolean flag for whether to print the documentation. default is False
    :param save_to_csv: A boolean flag for whether to save the matching results into a csv file. default is False.
    :param calculate_time: A boolean flag for whether to return the matching time for patient. default is False.
    If one wishes to save the results to csv files, a directory named 'Matching_Results_{searchId}' will be created in
    the working directory. If a directory by this name was already created, an error will be raised.
    Note: saving a pandas into a csv might take a couple of seconds
    :return: A dictionary that maps each patient to its matching results formatted as a pandas.DataFrame
    """
    if save_to_csv:
        os.makedirs(f"Matching_Results_{search_id}", exist_ok=True)

    if verbose:
        print_time("Start graph matching")

    g_m = DonorsMatching(match_graph, verbose=verbose)

    # create patients graph and find all candidates
    start_build_graph = time.time()
    subclasses_by_patient, classes_by_patient = g_m.create_patients_graph(imputation_filename)
    patients = list(g_m.patients.keys())

    if verbose:
        print_time("Created patients graph")
    end_build_graph = time.time()

    # the returned dictionary. {patient ID: pd.DataFrame(matches + features)}
    patients_results = {patient: None for patient in patients}

    if patients:
        avg_build_time = (end_build_graph - start_build_graph) / len(patients)
    else:
        avg_build_time = 0

    for patient in patients:
        # For each patient we search matches in the donor graph.
        # First we will look for perfect matches - only genotypes, then 9 matches - classes,
        # and then 7-8 matches - subclasses.
        if verbose:
            print_time(f"Searching matches for {patient}")

        start = time.time()

        subclasses = subclasses_by_patient[patient]
        classes = classes_by_patient[patient]
        results_df = search_in_levels(patient, g_m, donors_info, threshold, cutof, classes, subclasses)

        end = time.time()
        patient_time = end - start + avg_build_time

        if calculate_time:
            patients_results[patient] = (results_df, patient_time)
        else:
            patients_results[patient] = results_df

        if save_to_csv:
            # save results to csv
            results_df.to_csv(os.path.join(f"Matching_Results_{search_id}", f"Patient_{patient}.csv"), index=True,
                              float_format='%.2f')
            if verbose:
                print_time(f"Saved Matching results for {patient} in "
                           f"{os.path.join(f'Matching_Results_{search_id}', f'Patient_{patient}.csv')}")

    return patients_results


def matching(match_graph: MatchingGraph, grim_config_file="",
             save_imputation: Union[bool, str, PathLike] = False,
             donors_info: Union[Iterable[str], None] = None, search_id: int = 0,
             threshold: float = 0.1, cutof: int = 50,
             verbose: bool = False, save_to_csv: bool = False):
    """
    A function that performs the patients imputation with the matching.
    The imputation is performed with GRIM algorithm.

    :param match_graph: A Graph object from grma.match
    :param grim_config_file: A path to configuration file for grim imputation. default is grim default config.
    :param save_imputation: A flag for whether to save the imputation results. default is False.
    Accepts boolean/str/PathLike values - False will not save a file,
    True will save a file named 'imputation{searchId}.csv' in the working directory.
    str/PathLike will create and save the given file as long as the file's directory is valid.
    Otherwise, the file will not be saved. default is False.
    :param donors_info: An iterable of fields from the database to include in the results. default is None.
    :param search_id: An integer identification of the search. default is 0.
    :param threshold: Minimal score value for a valid match. default is 0.1.
    :param cutof: Maximum number of matches to return. default is 50.
    :param verbose: A boolean flag for whether to print the documentation. default is False
    :param save_to_csv: A boolean flag for whether to save the matching results into a csv file. default is False.
    If one wishes to save the results to csv files, a directory named 'Matching_Results_{searchId}' will be created in
    the working directory. If a directory by this name was already created, an error will be raised.
    Note: saving a pandas into a csv might take a couple of seconds
    :return: A dictionary that maps each patient to its matching results formatted as a pandas.DataFrame
    """
    if donors_info is None:
        donors_info = []

    # disable output from grim
    if not verbose:
        sys.stdout = open(os.devnull, 'w')

    imputation_path: str = run_grim(grim_config_file)

    # enable output if disabled
    if not verbose:
        sys.stdout = sys.__stdout__

    all_matches: Dict[int, pd.DataFrame] = find_matches(imputation_path, match_graph, search_id, donors_info,
                                                        threshold, cutof, verbose, save_to_csv)

    return all_matches
