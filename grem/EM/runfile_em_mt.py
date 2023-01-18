import argparse

# import cProfile
import json
import pathlib

import sys
import os
import numpy as np

# import threading
from multiprocessing import Process
import string, math
import copy

sys.path.insert(0, os.path.join(".."))


from grim import grim


# Profiler start
def run(
    plan_b,
    config_file="../conf/minimal-configuration.json",
    count_by_prob=None,
    pop=None,
    em_mr=False,
    iteration=1,
    num_subjects=None,
    project_dir_graph="",
    project_dir_in_file="",
):

    project_dir = ""

    # project_dir = ""# "../"
    output_dir = "output/"

    # Read configuration file and load properties
    with open(config_file) as f:
        json_conf = json.load(f)

    graph_files_path = json_conf.get("graph_files_path")
    if graph_files_path[-1] != "/":
        graph_files_path += "/"
    # output_dir = json_conf.get("imuptation_out_path", "output")
    # output_dir = "output/"
    # Read configuration file and load properties

    config = {
        "planb": json_conf.get("planb", True),
        "pops": json_conf.get("populations"),
        "priority": json_conf.get("priority"),
        "epsilon": json_conf.get("epsilon", 1e-3),
        "number_of_results": json_conf.get("number_of_results", 1000),
        "number_of_pop_results": json_conf.get("number_of_pop_results", 100),
        "output_MUUG": json_conf.get("output_MUUG", False),
        "output_haplotypes": json_conf.get("output_haplotypes", True),
        "node_file": project_dir_graph
        + graph_files_path
        + json_conf.get("node_csv_file"),
        "top_links_file": project_dir_graph
        + graph_files_path
        + json_conf.get("top_links_csv_file"),
        "edges_file": project_dir_graph
        + graph_files_path
        + json_conf.get("edges_csv_file"),
        "imputation_input_file": project_dir_in_file
        + json_conf.get("imputation_in_file"),
        "imputation_out_umug_freq_file": output_dir
        + json_conf.get("imputation_out_umug_freq_filename", "None"),
        "imputation_out_umug_pops_file": output_dir
        + json_conf.get("imputation_out_umug_pops_filename", "None"),
        "imputation_out_hap_freq_file": output_dir
        + json_conf.get("imputation_out_hap_freq_filename"),
        "imputation_out_hap_pops_file": output_dir
        + json_conf.get("imputation_out_hap_pops_filename", "None"),
        "imputation_out_miss_file": output_dir
        + json_conf.get("imputation_out_miss_filename", "miss.txt"),
        "imputation_out_problem_file": output_dir
        + json_conf.get("imputation_out_problem_filename", "problem.txt"),
        "factor_missing_data": json_conf.get("factor_missing_data", 0.01),
        "loci_map": json_conf.get(
            "loci_map", {"A": "A", "B": "B", "C": "C", "DQB1": "Q", "DRB1": "R"}
        ),
        "loci_map": json_conf.get(
            "loci_map", {"A": 1, "B": 3, "C": 2, "DQB1": 4, "DRB1": 5}
        ),
        "matrix_planb": json_conf.get(
            "Plan_B_Matrix",
            [
                [[1, 2, 3, 4, 5]],
                [[1, 2, 3], [4, 5]],
                [[1], [2, 3], [4, 5]],
                [[1, 2, 3], [4], [5]],
                [[1], [2, 3], [4], [5]],
                [[1], [2], [3], [4], [5]],
            ],
        ),
        "pops_count_file": project_dir + json_conf.get("pops_count_file", ""),
        "use_pops_count_file": json_conf.get("pops_count_file", False),
        "number_of_options_threshold": json_conf.get(
            "number_of_options_threshold", 100000
        ),
        "max_haplotypes_number_in_phase": json_conf.get(
            "max_haplotypes_number_in_phase", 100
        ),
        "bin_imputation_input_file": project_dir
        + json_conf.get("bin_imputation_in_file", "None"),
        "num_thread": json_conf.get("num_threads", 1),
        "nodes_for_plan_A": json_conf.get("Plan_A_Matrix", []),
        "save_mode": json_conf.get("save_space_mode", False),
        "UNK_priors": json_conf.get("UNK_priors", "MR"),
    }
    all_loci_set = set()
    for _, val in config["loci_map"].items():
        all_loci_set.add(str(val))
    config["full_loci"] = "".join(sorted(all_loci_set))
    config["imputation_out_hap_pops_file"] = (
        config["imputation_out_hap_pops_file"] + str(iteration) + ".txt"
    )

    if pop:
        config["pops"] = pop

    # Display the configurations we are using
    print(
        "****************************************************************************************************"
    )
    print("Performing imputation based on:")
    print("\tPopulation: {}".format(config["pops"]))
    print("\tPriority: {}".format(config["priority"]))
    print("\tEpsilon: {}".format(config["epsilon"]))
    print("\tPlan B: {}".format(config["planb"]))
    print("\tNumber of Results: {}".format(config["number_of_results"]))
    print("\tNumber of Population Results: {}".format(config["number_of_pop_results"]))
    print("\tNodes File: {}".format(config["node_file"]))
    print("\tTop Links File: {}".format(config["edges_file"]))
    print("\tInput File: {}".format(config["imputation_input_file"]))
    print("\tOutput UMUG Format: {}".format(config["output_MUUG"]))
    print(
        "\tOutput UMUG Freq Filename: {}".format(
            config["imputation_out_umug_freq_file"]
        )
    )
    print(
        "\tOutput UMUG Pops Filename: {}".format(
            config["imputation_out_umug_pops_file"]
        )
    )
    print("\tOutput Haplotype Format: {}".format(config["output_haplotypes"]))
    print(
        "\tOutput HAP Freq Filename: {}".format(config["imputation_out_hap_freq_file"])
    )
    print(
        "\tOutput HAP Pops Filename: {}".format(config["imputation_out_hap_pops_file"])
    )
    print("\tOutput Miss Filename: {}".format(config["imputation_out_miss_file"]))
    print("\tOutput Problem Filename: {}".format(config["imputation_out_problem_file"]))
    print("\tFactor Missing Data: {}".format(config["factor_missing_data"]))
    print("\tLoci Map: {}".format(config["loci_map"]))
    print("\tPlan B Matrix: {}".format(config["matrix_planb"]))
    print("\tPops Count File: {}".format(config["pops_count_file"]))
    print("\tUse Pops Count File: {}".format(config["use_pops_count_file"]))
    print(
        "\tNumber of Options Threshold: {}".format(
            config["number_of_options_threshold"]
        )
    )
    print(
        "\tMax Number of haplotypes in phase: {}".format(
            config["max_haplotypes_number_in_phase"]
        )
    )
    if config["nodes_for_plan_A"]:
        print("\tNodes in plan A: {}".format(config["nodes_for_plan_A"]))
    print("\tSave space mode: {}".format(config["save_mode"]))
    print(
        "****************************************************************************************************"
    )

    # Perform imputation
    graph = grim.graph_instance(config)
    imputation_list = []

    # Create output directory if it doesn't exist
    pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

    input_file = config["imputation_input_file"]
    in_dir = os.path.dirname(input_file)
    if in_dir == "":
        in_dir = "."
    if project_dir_in_file != "":
        in_dir = "splited_data_for_em"
        pathlib.Path(in_dir).mkdir(parents=False, exist_ok=True)
    in_file_basename = os.path.basename(input_file)

    if not num_subjects:
        num_subjects = os.popen("wc -l " + input_file).read()
        num_subjects = int(num_subjects.strip().split(" ")[0]) + 1
    print(num_subjects)

    split_cmd = (
        "split  -l "
        + str(int(math.ceil(num_subjects / config["num_thread"])))
        + " "
        + input_file
        + " "
        + in_dir
        + "/"
        + in_file_basename[0:2]
    )

    os.system(split_cmd)

    alpha = string.ascii_lowercase
    therads_list = list()
    config_list = []
    for i in range(config["num_thread"]):
        imputation = grim.impute_instance(config, graph, count_by_prob=count_by_prob)
        imputation_list.append(imputation)
        in_file = (
            in_dir
            + "/"
            + in_file_basename[0:2]
            + alpha[int(i / 26)]
            + alpha[int(i % 26)]
        )
        print(in_file)
        output_file = output_dir + "/" + os.path.basename(in_file) + "_out"
        config_list.append(copy.deepcopy(config))
        config_list[i]["imputation_input_file"] = in_file
        config_list[i]["imputation_out_hap_freq_file"] = output_file

    for i in range(config["num_thread"]):
        t = Process(
            target=imputation_list[i].impute_file,
            args=(
                config_list[i],
                plan_b,
                em_mr,
                True,
            ),
        )
        therads_list.append(t)
    for i in range(config["num_thread"]):
        therads_list[i].start()
    for i in range(config["num_thread"]):
        therads_list[i].join()
    for i in range(config["num_thread"]):
        therads_list[i].terminate()

    f_out = open(config["imputation_out_hap_freq_file"], "w")
    for i in range(config["num_thread"]):

        with open(config_list[i]["imputation_out_hap_freq_file"]) as f_t_out:
            for line in f_t_out:
                f_out.write(line)
            f_t_out.close()
            os.remove(config_list[i]["imputation_out_hap_freq_file"])

    # Profiler end
    # pr.disable()
    # pr.print_stats(sort="time")
