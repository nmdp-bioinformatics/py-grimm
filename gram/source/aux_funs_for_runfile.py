import os
import json

from .adding_children_data import add_child_data
from .parents_exist import associate_children_to_par_haps_while_parents_exist
from .als_update import Als


def get_files_path(out_files_path):
    out_GLstr = os.path.join(out_files_path, 'glstring_for_GRIMM.txt')
    out_binary = os.path.join(out_files_path, 'binary_for_GRIMM.txt')

    return out_GLstr, out_binary


def load_jsons():
    # cur_path = os.path.abspath("GR_code/GG_GRAMM/data")
    cur_path = os.path.dirname(os.path.realpath(__file__)) + '/data'

    with open(os.path.join(cur_path, 'low2high.txt')) as json1:
        low2high = json.load(json1)

    with open(os.path.join(cur_path, 'ser_dict_antigen2group.json')) as json2:
        antigen2group = json.load(json2)

    with open(os.path.join(cur_path, 'ser_dict_group2antigen.json')) as json3:
        group2antigen = json.load(json3)

    return low2high, antigen2group, group2antigen


def check_num_parents(fam_dict):
    """
    check how many parents exist in the family
    :param fam_dict: family dict
    :return: number of parents
    """
    fam_members = list(fam_dict.keys())
    if 'F' in fam_members and 'M' in fam_members:
        return 2
    if 'F' in fam_members or 'M' in fam_members:
        return 1
    return 0


def convert_data_to_Als(fam_dict):
    """
    convert, for a family, the alleles_names data format: from a list to an Als
    (like list, just adjusted to alleles_names. see 'Als' class documentation)
    :param fam_dict: family dict
    """
    for fam_member in fam_dict:  # F, M, 1 ...
        for allele_name in fam_dict[fam_member]:  # A, B ...
            al1 = fam_dict[fam_member][allele_name][0]  # first allele
            al2 = fam_dict[fam_member][allele_name][1]  # second allele
            new_format = Als()  # new object of Als
            new_format.extend([al1, al2])  # add alleles_names data

            fam_dict[fam_member].update({allele_name: new_format})  # update the data in the dict to be in Als format


def remove_duplicate_children(family, alleles_names):
    """
    if there are two children with identical alleles_names data, remove one
    :param family: family dict
    :param alleles_names: alleles_names
    """

    def duplicate(child1, child2):
        for al_name in alleles_names:
            data1 = family[child1][al_name]
            data2 = family[child2][al_name]
            if data1 != data2:
                return False
        return True

    tested_children = []  # insert to this list the children we checked
    to_remove = []  # children to remove. (don't del them into the loop because it crash,
    # cause can't change dict while we run on it)
    for child in family:
        if child not in ['F', 'M']:
            if not tested_children:  # empty list
                tested_children.append(child)
            else:
                for other_child in tested_children:
                    if other_child != child:  # check it's not the same _child_
                        are_identical = duplicate(child, other_child)
                        if are_identical:
                            to_remove.append(child)
                            break
                        else:
                            tested_children.append(child)

    for dup_child in to_remove:
        del family[dup_child]


def insert_parents_data_to_their_haps(family, hapF, hapM, par_num):
    """
    when parents exist (al least one), we insert their data to their haplotypes
    (explanation about the insertion manner in documentation of 'insert_parents_data')
    :param family: family dict
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param par_num: how many parents (0/1/2)
    """
    if 'F' in family:
        hapF.insert_parents_data(family, 'F', par_num)
    if 'M' in family:
        hapM.insert_parents_data(family, 'M', par_num)


def do_planB(hapF, hapM, child, child_idx, children_num, alleles_names, aux_tools, count_fam, errors_in_families):
    """
    trying to associate and add data again from children that the associating before wasn't complete
    - we send to the function that the parents exist ('associate_children...parents_exist'), because we now after
      loop over the other children, so probably will be data in the parents haplotypes
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param child: _child_ data (dict)
    :param child_idx: _child_ index
    :param children_num: number of children in family
    :param alleles_names: alleles_names names
    :param aux_tools: dict with auxiliary tools
    :param count_fam: family index
    :param errors_in_families: dict with data about error in the families
    :return: if the adding is succeed or not
    """
    success_associate, associating = \
        associate_children_to_par_haps_while_parents_exist(hapF, hapM, child, alleles_names)
    success_adding = add_child_data(hapF, hapM, child, child_idx, children_num, associating, alleles_names, aux_tools,
                                    count_fam, errors_in_families)
    return success_adding


def write_binaries_to_file(out_binary, aux_tools):
    """
    write binary list of each pearson to a file, for GRIMM
    :param out_binary: out file path
    :param aux_tools: dict which contains the binaries dict
    """
    bin_dict = aux_tools['binary_dict']
    if len(bin_dict) > 0:
        with open(out_binary, 'w') as file:
            json.dump(bin_dict, file)






