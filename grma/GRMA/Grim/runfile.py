import os
import warnings
from typing import Union

from grim.grim import impute_instance
from grim.imputation.networkx_graph import Graph
from GRMA.Grim.configuration import CONFIG


def run_GRIMM(inputFile: Union[str, os.PathLike], graph: Graph,
              searchId: int, output_dir: Union[str, os.PathLike], save_imputation: Union[bool, str, os.PathLike]):
    """
    A function for running the grim imputation algorithm.

    :param inputFile: donors filepath
    :param graph: A Graph object from grim.imputation.networkx_graph
    :param searchId: An integer identification of the search.
    :param output_dir: A PathLike for the directory of the output files of grim.
    :param save_imputation: A flag for whether to save the imputation results.
    Accepts boolean/str/PathLike values - False will not save a file,
    True will save a file named 'imputation{searchId}.csv' in the working directory.
    str/PathLike will create and save the given file as long as the file's directory is valid.
    Otherwise, the file will not be saved. default is False.
    :return: A PathLike for the created imputation file.
    """

    # validate the input file location
    if not (os.path.exists(inputFile) and os.path.isfile(inputFile)):
        raise FileNotFoundError(f"inputFile: {inputFile} for grim cannot be fount")

    # validata the given input of imputation_path
    if isinstance(save_imputation, bool):
        imputation_path = os.path.join(output_dir, f"imputations{searchId}.csv") if not save_imputation \
            else f"imputations{searchId}.csv"
    else:
        splited = os.path.split(save_imputation)
        if os.path.exists(splited[0]) and os.path.isdir(splited[1]):
            imputation_path = save_imputation
        else:
            warnings.warn(f"Cannot find the given directory of the imputation file {save_imputation}. "
                          f"DO NOT SAVE THE IMPUTATION FILE.", ResourceWarning)
            imputation_path = os.path.join(output_dir, f"imputations{searchId}.csv")

    # set the rest of the configuration
    config = CONFIG
    config["imputation_input_file"] = inputFile
    config["imputation_out_umug_freq_file"] = imputation_path
    config["imputation_out_umug_pops_file"] = os.path.join(output_dir, f"out.umug.pops")
    config["imputation_out_hap_freq_file"] = os.path.join(output_dir, f"out.freqs")
    config["imputation_out_hap_pops_file"] = os.path.join(output_dir, f"out.pops")
    config["imputation_out_miss_file"] = os.path.join(output_dir, f"out.miss")
    config["imputation_out_problem_file"] = os.path.join(output_dir, f"out.problem")

    all_loci_set = set()
    for _, val in config["loci_map"].items():
        all_loci_set.add(str(val))

    config["full_loci"] = ''.join(sorted(all_loci_set))

    # impute the input file.
    imputation = impute_instance(config, graph)

    # Write out the results from imputation
    imputation.impute_file(config)

    return imputation_path
