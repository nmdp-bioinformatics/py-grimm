import os
import shutil
import sys
from os import PathLike
from typing import Iterable, Union, Dict

import pandas as pd
from grim.imputation.networkx_graph import Graph as GrimGraph

from GRMA.Grim import run_GRIMM
from GRMA.Match import Graph as MatchingGraph
from GRMA.Match.DonorsMatching import DonorsMatching, _init_results_df
from GRMA.Match.GraphWrapper import Graph
from GRMA.Utilities.utils import print_time


def find_matches(imputation_filename: Union[str, PathLike], graph: Graph,
                 searchId: int, donors_info: Iterable[str],
                 threshold: float = 0.1, cutof: int = 50,
                 verbose: bool = False, save_to_csv: bool = False) -> Dict[int, pd.DataFrame]:
    """
    The main function responsible for performing the matching.
    Note: for each patient, if a donor has been found as a
    match in an early stage (0, 1, or 2 mm), he will not be searched as a match for the further mismatches.

    :param imputation_filename: Path to the output file of the imputation made by grim.
    :param graph: A Graph object from GRMA.Match
    :param searchId: An integer identification of the search. default is 0.
    :param donors_info: An iterable of fields from the database to include in the results. default is None.
    :param threshold: Minimal score value for a valid match. default is 0.1.
    :param cutof: Maximum number of matches to return. default is 50.
    :param verbose: A boolean flag for whether to print the documentation. default is False
    :param save_to_csv: A boolean flag for whether to save the matching results into a csv file. default is False.
    If one wishes to save the results to csv files, a directory named 'Matching_Results_{searchId}' will be created in
    the working directory. If a directory by this name was already created, an error will be raised.
    Note: saving a pandas into a csv might take a couple of seconds
    :return: A dictionary that maps each patient to its matching results formatted as a pandas.DataFrame
    """
    if save_to_csv:
        os.mkdir(f"Matching_Results_{searchId}")

    if verbose:
        print_time("Start graph matching")

    g_m = DonorsMatching(graph)

    # create patients graph and find all candidates
    patients = g_m.create_patients_graph_and_find_candidates(imputation_filename, verbose)
    if verbose:
        print_time("Created patients graph")

    # the returned dictionary. {patient ID: pd.DataFrame(matches + features)}
    patients_results = {patient: None for patient in patients}

    # loop over all patients
    for j, pat in enumerate(patients.keys()):
        if verbose:
            print_time(f"Searching matches for {pat}")

        matched = set()  # set of donors ID that have already matched for this patient
        count_matches = 0
        results_df = _init_results_df(donors_info)  # initialize the df according to the given fields

        # loop over possible mismatches: 0, 1, 2, 3.
        for mismatches in range(4):
            matched, count, results_df = g_m.score_matches(mismatches, results_df, donors_info,
                                                           pat, threshold, cutof, matched)
            count_matches += count
            if verbose:
                print_time(f"({mismatches} MMs) Found {count} matches")

        patients_results[pat] = results_df

        if save_to_csv:
            # save results to csv
            results_df.to_csv(os.path.join(f"Matching_Results_{searchId}", f"Patient_{pat}.csv"), index=True,
                              float_format='%.2f')
            if verbose:
                print_time(f"Saved Matching results for {pat} in "
                           f"{os.path.join(f'Matching_Results_{searchId}', f'Patient_{pat}.csv')}")

    return patients_results


def matching(filepath: Union[str, PathLike], grim_graph: GrimGraph, match_graph: MatchingGraph,
             save_imputation: Union[bool, str, PathLike] = False,
             donors_info: Union[Iterable[str], None] = None, searchId: int = 0,
             threshold: float = 0.1, cutof: int = 50,
             verbose: bool = False, save_to_csv: bool = False) -> Dict[int, pd.DataFrame]:
    """
    A function that performs the patients imputation with the matching.
    The imputation is performed with GRIM algorithm.

    :param filepath: path to patients file.
    :param grim_graph: A Graph object from grim.imputation.networkx_graph
    :param match_graph: A Graph object from GRMA.Match
    :param save_imputation: A flag for whether to save the imputation results. default is False.
    Accepts boolean/str/PathLike values - False will not save a file,
    True will save a file named 'imputation{searchId}.csv' in the working directory.
    str/PathLike will create and save the given file as long as the file's directory is valid.
    Otherwise, the file will not be saved. default is False.
    :param donors_info: An iterable of fields from the database to include in the results. default is None.
    :param searchId: An integer identification of the search. default is 0.
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

    directory = f"MatchingFiles{searchId}"
    try:
        # create a directory for the output of grim
        os.makedirs(directory, exist_ok=True)

        # disable output from grim
        if not verbose:
            sys.stdout = open(os.devnull, 'w')

        imputation_filename: str = run_GRIMM(filepath, grim_graph, searchId, directory, save_imputation)

        # enable output if disabled
        if not verbose:
            sys.stdout = sys.__stdout__

        all_matches: Dict[int, pd.DataFrame] = find_matches(imputation_filename, match_graph, searchId, donors_info,
                                                            threshold, cutof, verbose, save_to_csv)
    finally:
        # delete the directory created for the output of grim
        if os.path.exists(directory) and os.path.isdir(directory):
            shutil.rmtree(directory)

    return all_matches
