import os
import json
from datetime import datetime
from shutil import move
from gramm.source.runfile_update import run_GRAMM
from gramm.source.create_results_updated import run_Post_GRIMM
from gramm.source.visualization_updated import visualize

from grim import grim


def gramm(config_file=None):
    if config_file is None:
        config_file = os.path.dirname(os.path.realpath(__file__)) + '/static/config_gramm.json'

    with open(config_file) as f:
        config = json.load(f)

    user_file = os.path.dirname(os.path.realpath(__file__)) + config.get('user_file_path', '/static/example_file1.csv')
    config_file_grimm = os.path.dirname(os.path.realpath(__file__)) + config.get('config_file_grimm', '/static/config_grim.json')
    output_dir = config.get('output_dir', 'user_files')
    alleles_names = config.get('alleles_names', ['A', 'B', 'C', 'DRB1', 'DQB1'])
    is_serology_data = config.get('is_serology_data', False)
    race_dict = config.get("races", {})

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    run_all(user_file, alleles_names, output_dir, is_serology_data, race_dict, config_file_grimm)
    organize_files(output_dir)


def run_all(input2GRAMM, alleles_names, files_address, is_serology, race_dict, config_grimm, open_ambiguity_sim=False):
    """
    run GRAMM, GRIMM, Post-GRAMM
    :param input2GRAMM: file with input from user
    :param alleles_names: ['A', 'B', 'C', 'DRB1', 'DQB1']
    :param files_address: path to directory (that stores the output files)
    :param is_serology: if data is serology, True. otherwise, False
    :param race_dict: dict that contains the races of the family that the user inserted
    :param config_grimm: config file for grim
    :param open_ambiguity_sim: in some simulations and Israeli-Australian file data, it keeps the data of the ambiguity
        of the families, that have been removed before, to use it in GRAMM
    :return: path to results file, path to errors file, flag that=True if the process did not find any valid results
    """

    errors_in_families, aux_tools = run_GRAMM(input2GRAMM, files_address, alleles_names, race_dict,
                                              open_ambiguity_sim, is_serology)

    grim.graph_freqs(conf_file=config_grimm)
    grim.impute(conf_file=config_grimm)

    input2post = f'{files_address}/output_grim/don.hap.freqs'
    results, run_again = run_Post_GRIMM(input2post, input2GRAMM, alleles_names, files_address,
                                        aux_tools, errors_in_families)
    visualize(results, aux_tools['parent_has_empty_hap'], files_address)

    return results, files_address + '/errors.txt', run_again


def organize_files(output_dir):
    curr_time = datetime.now().strftime("%Y%m%d_%H%M")
    output_from_gramm = f'{output_dir}/output_gramm_{curr_time}'
    os.mkdir(output_from_gramm)

    for filename in os.listdir(output_dir):
        file_path = os.path.join(output_dir, filename)
        if os.path.isfile(file_path):
            if filename in ['glstring_for_GRIMM.txt', 'binary_for_GRIMM.txt']:
                os.remove(file_path)
            else:
                move(file_path, os.path.join(output_from_gramm, filename))