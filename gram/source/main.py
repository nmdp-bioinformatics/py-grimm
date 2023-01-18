import os
import sys
import json
from datetime import datetime
from pathlib import Path
from .runfile_update import run_GRAMM
from .create_results_updated import run_Post_GRIMM
from .visualization_updated import visualize

from grim import grim

# sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)).replace("/source", ""))


def gram(config_file=None):
    if config_file is None:
        config_file = Path(__file__).parents[1] / Path('conf/gram_config.json')

    with open(config_file) as f:
        config = json.load(f)

    # get gram parameters
    input_file_path = config.get('input_file_path', 'gram/source/data/input_file_example.csv')
    output_dir = config.get('output_dir', 'output_gram')
    alleles_names = config.get('alleles_names', ['A', 'B', 'C', 'DRB1', 'DQB1'])
    is_serology_data = config.get('is_serology_data', False)
    race_dict = config.get('races', {})

    # get grim parameters
    config_grim_path = config.get('config_grim_path', '')
    build_grim_graph = config.get('build_grim_graph', True)

    output_dir = create_output_dir(output_dir)
    run_all(input_file_path, alleles_names, output_dir, is_serology_data, race_dict, config_grim_path, build_grim_graph)
    remove_unnecessary_files(output_dir)


def create_output_dir(output_dir):
    curr_time = datetime.now().strftime("%Y%m_%d%H%M")
    output_dir = f'{output_dir}_{curr_time}'
    os.mkdir(output_dir)
    return output_dir


def run_all(input_file_path, alleles_names, output_dir, is_serology, race_dict, grim_config_path, build_grim_graph,
            open_ambiguity_sim=False):
    """
    run 3 parts of gram (gram, grim, post-gram)
    :param input_file_path: file with input from user
    :param alleles_names: ['A', 'B', 'C', 'DRB1', 'DQB1']
    :param output_dir: output directory path
    :param is_serology: whether the input data is serology
    :param race_dict: a dictionary that contains the races of the families that the user inserted
    :param grim_config_path: path the grim configuration file
    :param build_grim_graph: flag - whether run the function 'grim.graph_freqs'
           (that build the frequencies graph for grimm). its required in the first running only.
    :param open_ambiguity_sim: relevalnt in simulation files
    """
    errors_in_families, aux_tools = run_GRAMM(input_file_path, output_dir, alleles_names, race_dict,
                                              open_ambiguity_sim, is_serology)
    if build_grim_graph:
        grim.graph_freqs(conf_file=grim_config_path)
    grim_conf_adjusted, out_imputation_file = process_grim_config_file(grim_config_path, output_dir)
    grim.impute(conf_file=grim_conf_adjusted)

    results, run_again = run_Post_GRIMM(out_imputation_file, input_file_path, alleles_names, output_dir,
                                        aux_tools, errors_in_families)
    visualize(results, aux_tools['parent_has_empty_hap'], output_dir)


def process_grim_config_file(config_grim_path, output_dir):
    """
    read grim config file, adjust it for gram imputation and write a adjusted config file to 'conf' directory.
    :param config_grim_path: path to grim config file
    :param output_dir: path to gram output directory
    :return: a path to the adjusted grim config file and a path to grim imputation output
    """
    if config_grim_path == '':
        config_grim_path = Path(__file__).parents[2] / Path('conf/minimal-configuration.json')
        if not os.path.exists(config_grim_path):  # in case we do not found the default config of grim, use a local copy
            config_grim_path = Path(__file__).parents[1] / Path('conf/grim_conf_copy.json')

    with open(config_grim_path) as config_grim:
        grim_dict = json.load(config_grim)

    grim_dict['imputation_in_file'] = os.path.join(output_dir, 'glstring_for_GRIMM.txt')
    grim_dict['bin_imputation_input_file'] = os.path.join(output_dir, 'binary_for_GRIMM.txt')

    out_imputation_file = os.path.join(grim_dict['imuptation_out_path'], grim_dict['imputation_out_hap_freq_filename'])

    grim_conf_adjusted = Path(__file__).parents[1] / Path('conf/grim_conf_adjusted.json')
    with open(grim_conf_adjusted, 'w') as new_grim_config:
        json.dump(grim_dict, new_grim_config)

    return grim_conf_adjusted, out_imputation_file


def remove_unnecessary_files(output_dir):
    for filename in os.listdir(output_dir):
        file_path = os.path.join(output_dir, filename)

        if os.path.isfile(file_path):
            if filename in ['glstring_for_GRIMM.txt', 'binary_for_GRIMM.txt']:
                os.remove(file_path)
