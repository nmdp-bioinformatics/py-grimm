import json
import os
import sys

sys.path.insert(0, os.path.join(".."))

from . import em_haplotypes_algo
from . import em_mr_algo
from graph_generation import generate_neo4j_multi_hpf
from . import runfile_em_mt
import numpy as np
import pathlib

# import shutil


def run_em_def(
    conf_file="../conf/minimal-em-configuration.json",
    sr_pop_name="all",
    project_dir_in_file="",
):

    project_dir = "../"
    output_dir = "output/"
    # photos_dir = 'photos/'
    graph_generation_dir = "../imputation/graph_generation/"
    # Read configuration file and load properties
    with open(conf_file) as f:
        json_conf = json.load(f)
    graph_files_path = json_conf.get("graph_files_path")
    config = {
        "imputation_input_file": project_dir_in_file
        + json_conf.get(
            "imputation_in_file"
        ),  # "/".join(json_conf.get("imputation_in_file").strip("/").split('/')[1:]) ,
        "freq_file": json_conf.get("freq_file"),
        "imputation_out_hap_freq_file": output_dir
        + json_conf.get("imputation_out_hap_freq_filename"),
        "pops": json_conf.get("populations"),
        "loci_map": json_conf.get("loci_map"),
        "cutoff_init": json_conf.get("init_cutoff", 100),
        "log_file_name": json_conf.get("logLikelihood_file", "log_likelihood.txt"),
        "memory_max": json_conf.get("memory_max_size", 6000000),
        "memory_min": json_conf.get("memory_min_size", 1000000),
        "max_iteration": json_conf.get("max_iterations", 50),
        "just_SR": json_conf.get("run_just_SR_EM", False),
        "node_file": graph_files_path + json_conf.get("node_csv_file"),
        "top_links_file": graph_files_path + json_conf.get("top_links_csv_file"),
        "edges_file": graph_files_path + json_conf.get("edges_csv_file"),
        "info_nodes": graph_files_path + json_conf.get("info_node_csv_file"),
    }

    file_lo = open(config["log_file_name"], "w")

    pop = sr_pop_name

    name_f_freq = os.path.basename(config.get("freq_file"))

    # Create output directory if it doesn't exist
    pathlib.Path("/".join(config.get("freq_file").split("/")[:-1])).mkdir(
        parents=False, exist_ok=True
    )
    pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

    ### em without races
    logL = -100
    algo = em_haplotypes_algo.algo(config, pop)

    eps, numH, num_saples = algo.create_guess()

    iteration = 0
    not_converge = True
    planB = False
    while not_converge and iteration < config["max_iteration"]:
        print("iter = " + str(iteration))

        # file_for_save = output_dir + ls + '_SR_res_' + str(iteration) + '.csv'
        # shutil.copy2(config["freq_file"], file_for_save)

        iteration += 1

        generate_neo4j_multi_hpf.generate_graph(conf_file, [pop], em=True)

        runfile_em_mt.run(
            plan_b=planB,
            config_file=conf_file,
            pop=[pop],
            count_by_prob=np.array([1]),
            num_subjects=num_saples,
            project_dir_in_file=project_dir_in_file,
        )
        last_logl = logL
        logL, numH = algo.calc_new_prob()
        logL /= num_saples
        file_lo.write(str(iteration) + " " + str(logL / num_saples) + "\n")
        planB = True
        if logL - last_logl <= 0.01 and iteration > 2:
            not_converge = False

    print(logL)

    """file_for_kl_mr = output_dir +  'SR_res' + ls + '.csv'
    count_file_for_kl_mr = sum(1 for line in open(config["freq_file"]))
    shutil.copy2(config["freq_file"], file_for_kl_mr)"""

    ### em - mr
    if not config["just_SR"]:
        eps = 0
        algo = em_mr_algo.algo_mr(config, eps)
        eps = algo.create_eps()
        algo.create_guess()

        not_converge = True
        planB = True
        iteration = 0
        pops_len = len(config["pops"])
        count_by_prob = np.ones(pops_len)

        while not_converge and iteration < config["max_iteration"]:
            print("iter = " + str(iteration))

            iteration += 1
            """dict_hf = {}
            with open(config["freq_file"]) as hpf_file:
                for line in hpf_file:
                    hap, race, prob = line.strip().split(',')
                    if not hap in dict_hf:
                        dict_hf[hap] = np.zeros(pops_len)
                    dict_hf[hap][config["pops"].index(race)] = float(prob)"""

            generate_neo4j_multi_hpf.generate_graph(conf_file, em=True)
            runfile_em_mt.run(
                planB,
                conf_file,
                count_by_prob=count_by_prob,
                em_mr=True,
                iteration=iteration,
                num_subjects=num_saples,
                project_dir_in_file=project_dir_in_file,
            )
            last_logl = logL
            count_by_prob, logL = algo.calc_new_prob()
            file_lo.write(str(iteration) + " " + str(logL) + "\n")
            planB = True
            # not_converge = algo.check_converges(dict_hf)
            if logL - last_logl <= 0.01:
                not_converge = False
        print("Log Likelihood: " + str(logL))

    os.remove(config["node_file"])
    os.remove(config["top_links_file"])
    os.remove(config["edges_file"])
    os.remove(config["info_nodes"])

    file_lo.write("loglikelihood " + str(logL) + "\n")
