import os
import csv
import json
import itertools

from .read_families_from_file import get_families
from .general_aux_funs import convert_to_serology
from .als_update import Als


def run_Post_GRIMM(input_from_grimm, orig_input, alleles_names, output_path, aux_tools, errors_in_families):
    """
    the purpose of "POST-GRIMM" is comparing the gl-string results of parents we get from GRIMM, with the data about
    children.
    each result (of pair haplotypes from father and pair from mother) that not consistent with the children - reject
    :param input_from_grimm: path to results file from GRIMM
    :param orig_input: the input from the user (after format adjustment)
    :param alleles_names: alleles names
    :param output_path: the path to output files (results and errors)
    :param aux_tools: dict with auxiliary tools
    :param errors_in_families: dict that contains the errors of families
    :return:
    """
    min_freq = 1e-15  # hyper-parameter for limit the frequency of results (only if there is already valid results)

    res_file, writer = open_res_file(output_path)

    input_from_grimm = open(input_from_grimm).readlines()
    grimm_dict = {}

    # create dict with the data from GRIMM
    for i in range(len(input_from_grimm)):
        input_from_grimm[i] = input_from_grimm[i].strip().split(",")
        if input_from_grimm[i][0] in grimm_dict:
            grimm_dict[input_from_grimm[i][0]].append(input_from_grimm[i][1:])
        else:
            grimm_dict[input_from_grimm[i][0]] = [input_from_grimm[i][1:]]

    families_dict = get_families(orig_input, aux_tools, dont_open_ser=True)  # get families from the original file (data from user)

    if aux_tools['is_serology']:  # serology data. open children to serology options
        antigen2group = aux_tools['antigen2group']
        group2antigen = aux_tools['group2antigen']

        for key_fam, family_dict in families_dict.items():
            for idx_member, member in family_dict.items():
                # if idx_member not in ['F', 'M']:  # open to serology only children data
                    for allele_name, alleles_values in member.items():
                        for idx_allele, single_allele in enumerate(alleles_values):
                            # we first convert to group name and then open all the options to the group name
                            families_dict[key_fam][idx_member][allele_name][idx_allele] = convert_to_serology(antigen2group,
                                                                                                 allele_name,
                                                                                                 single_allele)
                    open_serology_options_children(member, alleles_names, group2antigen)

    valid, invalid = 0, 0
    valid_list, invalid_list = [], []

    # compare between original data about families and the data from GRIMM
    # by checking if children (from input) are consistent with parents haplotypes (from GRIMM).
    for key in families_dict:

        if key + '.0' in grimm_dict.keys() and key + '.1' in grimm_dict.keys():
            p1_lst = grimm_dict[key + '.0']  # options to glstring of parent 1
            p2_lst = grimm_dict[key + '.1']  # options to glstring of parent 2
            # convert freq to float, then sort the list (for starting with the max freq)
            for p in p1_lst:
                p[2] = float(p[2])
            for p in p2_lst:
                p[2] = float(p[2])
            p1_lst = sorted(p1_lst, key=lambda x: (x[2]))
            p2_lst = sorted(p2_lst, key=lambda x: (x[2]))
            p1_lst = p1_lst[::-1]
            p2_lst = p2_lst[::-1]

            max_res2fam = 0
            for pair1 in p1_lst:  # we go over by nested loop on all the pairs combination between father and mother
                if max_res2fam == 100:  # max results for a family: 100
                    break

                # check consistency between parent (F) and all children
                if not consistent_haps_par_and_geno_child(pair1, families_dict[key], alleles_names):
                    # incons_haps += 1
                    continue

                for pair2 in p2_lst:
                    if max_res2fam == 100:
                        break  # max results of family is 100

                    # check consistency between parent (M) and all children
                    if not consistent_haps_par_and_geno_child(pair1, families_dict[key], alleles_names):
                        continue

                    # if fam key already appears in results and the freqs is too low, it won't be added to results:
                    if key not in valid_list or float(pair1[2]) * float(pair2[2]) > min_freq:
                        if valid_family(key, families_dict[key], pair1, pair2, alleles_names, aux_tools,
                                        errors_in_families, writer):
                            max_res2fam += 1
                            if key not in valid_list:
                                valid += 1
                                valid_list.append(key)
            if key not in valid_list:  # no combination of parents haplotypes was consistent with children
                invalid += 1
                invalid_list.append(key)
                errors_in_families[key] = ['All', '9']

    errors_summary_and_writing_to_file(output_path, errors_in_families, len(families_dict))

    res_file.close()

    # run_again is a flag that give sign if there is 1 family and it failed_count because inconsistent.
    # in this case, we run again the process, with 1000 results in grimm (instead of 10).
    # it can lead to results (when the alleles are rare)
    run_again = True if len(families_dict) == 1 and len(invalid_list) == 1 else False
    return output_path + '/results.csv', run_again


def open_res_file(output_path):
    res_file = open(output_path + '/results.csv', 'w', newline='')
    writer = csv.writer(res_file)
    row = ['FAMCODE', 'HAP_F1', 'HAP_F2', 'HAP_M1', 'HAP_M2', 'HAPS_INHERITANCE', 'PROBABILITY']
    writer.writerow(row)
    return res_file, writer


def open_serology_options_children(child, alleles_names, s_dict):
    """
    convert name of serology group to all the antigens in this group
    :param child: _child_ data
    :param alleles_names: alleles names
    :param s_dict: serologu dict
    """
    for i, alleles in enumerate(child):
        for al in child[alleles]:
            group = alleles_names[i] + '*' + al
            # check if the value (group) is in dict (low res only)
            # if it exists, open to all the serology options and remove the original value
            # (it's only group name, not an antigen)
            # example: A*10 appears, so add A*25, A*25, A*34, A*66 (only the digits) and remove A*10
            if ":" not in al and group in s_dict:
                for ser_value in s_dict[group]:
                    child[alleles].append(ser_value.split('*')[1])


def consistent_haps_par_and_geno_child(pair, family, alleles_names):
    """
    compare single parent to all the children and check consistency
    (it's a necessary condition for consistency with two parents)
    :param pair: pair of rwo haplotypes
    :param family: family dict
    :param alleles_names: alleles names
    :return: True if consistency, False otherwise
    """
    hap1 = gl_string_to_dict(pair[0])
    hap2 = gl_string_to_dict(pair[1])

    for idx_member, member in family.items():
        if idx_member not in ['F', 'M']:
            for allele_name in alleles_names:
                if not common_allele([hap1[allele_name], hap2[allele_name]], member[allele_name]):
                    return False
    return True


def gl_string_to_dict(gl_string):
    """
    from 'A*02:01~B*40:01..' to {A: 02:02, B:40:01 ..}
    """
    return {allele.split('*')[0]: allele.split('*')[1] for allele in gl_string.split('~')}


def common_allele(als_list1, als_list2):
    """
    check if exist common alleles between two list
    [01:02, 03:05] ; [01:03, 04:05] -- > return False because there is no common allele
    [01:02, 03:05] ; [01:02, 04:05] -- > return True because 01:02 is common
    :param als_list1: first list of alleles
    :param als_list2: second list pf alleles
    :return: True if exist, False otherwise
    """
    for al1 in als_list1:
        for al2 in als_list2:
            if is_equal(al1, al2):
                return True
    return False


def sort_hap(hap):  # organize in this order: A,B,C,DRB1,DQB1
    hap = hap.split('~')
    hap[-1], hap[-2] = hap[-2], hap[-1]
    hap = '~'.join(hap)
    return hap


def valid_family(key, family, p_1, p_2, alleles_names, aux_tools, errors_in_families, writer):
    """
    check if a given combination of parents haplotypes is consistent with all children
    :param key: family key (index)
    :param family: family dict
    :param p_1: pair haplotypes of parent 1
    :param p_2: pair haplotypes of parent 2
    :param alleles_names: alleles names
    :param aux_tools: dict with auxiliary tools
    :param errors_in_families: dict that contains families with errors
    :param writer: writer to results file
    :return: True if consistency, False otherwise
    """
    hap_1, hap_2, hap_3, hap_4 = sort_hap(p_1[0]), sort_hap(p_1[1]), sort_hap(p_2[0]), sort_hap(p_2[1])
    freq_1, freq_2 = float(p_1[2]), float(p_2[2])
    valid, not_valid = [], []
    is_serology = aux_tools['is_serology']

    for idx_member, member in family.items():
        if idx_member in ['F', 'M']:
            if validate(hap_1, hap_2, member, is_serology) or validate(hap_3, hap_4, member, is_serology):
                valid.append(idx_member)
            else:
                not_valid.append(idx_member)
        else:  # a _child_
            if validate(hap_1, hap_3, member, is_serology):
                valid.append(''.join(['C', idx_member, '=F1~M1']))
            elif validate(hap_1, hap_4, member, is_serology):
                valid.append(''.join(['C', idx_member, '=F1~M2']))
            elif validate(hap_2, hap_3, member, is_serology):
                valid.append(''.join(['C', idx_member, '=F2~M1']))
            elif validate(hap_2, hap_4, member, is_serology):
                valid.append(''.join(['C', idx_member, '=F2~M2']))
            else:
                not_valid.append(idx_member)

    # second condition for allow one problematic child. (but only if there are at least 2 children)
    if len(not_valid) == 0 or (len(not_valid) == 1 and not_valid[0] not in ['F', 'M'] and any(['~' in v in valid for v in valid])):
        match_child2haps = [member for member in valid if '~' in member]  # example: [1=F2~M2, 2=F1~M2 ...]

        if key in errors_in_families and errors_in_families[key][0] != 'All':
            idx_problematic_child = errors_in_families[key][0]
            # _child_ with error (but the rest of the family is valid)
            match_child2haps.append('C' + idx_problematic_child + '=XX~XX')

        match_child2haps = ";".join(match_child2haps)

        # sign to parent that had only one haplotype with data, so we duplicate it before grimm
        # (because grimm can't get single haplotype), and now we remove the duplicate and sign as "unknown"
        parent_has_empty_hap = aux_tools['parent_has_empty_hap']
        if key in parent_has_empty_hap:
            if parent_has_empty_hap[key] == 'F':
                hap_2 = 'Unknown'
            elif parent_has_empty_hap[key] == 'M':
                hap_4 = 'Unknown'

        row = [key, hap_1, hap_2, hap_3, hap_4, match_child2haps, freq_1 * freq_2]
        writer.writerow(row)
        return True

    return False


def validate(hap_1, hap_2, member, is_serology):
    """
    compare two haplotype to person (2 from one parent if compare to parent, and 1 from each parent if compare to _child_)
    :param hap_1: first haplotype
    :param hap_2: second haplotype
    :param member: family member
    :param is_serology: flag, if serology, the validate is checked in another way
    :return: True if consistency, False otherwise
    """
    hap_1 = gl_string_to_dict(hap_1)
    hap_2 = gl_string_to_dict(hap_2)

    if not is_serology:
        for allele_name, allele_values in member.items():
            val_member1, val_member2, val_hap1, val_hap2 = \
                allele_values[0], allele_values[1], hap_1[allele_name], hap_2[allele_name]
            success_option1 = is_equal(val_member1, val_hap1) and is_equal(val_member2, val_hap2)
            success_option2 = is_equal(val_member1, val_hap2) and is_equal(val_member2, val_hap1)
            if not(success_option1 or success_option2):
                return False
        return True

    else:  # serology data. could be more than 2 options in 'member[allele_name]'
        pairs_consistent = [False] * len(member.keys())
        for idx_allele, (allele_name, allele_values) in enumerate(member.items()):
            val_hap1, val_hap2 = hap_1[allele_name], hap_2[allele_name]
            val_member = Als()
            val_member.extend(member[allele_name])

            # the first two conditions are for [] and ["", ""]
            if not val_member or not any(val_member) or val_hap1 in val_member and val_hap2 in val_member:
                # avoid cases like '01' in val_member once, and val_hap1='01:01', val_hap2='01:02' (so it's invalid, but previous condition is true)
                if not (val_member.index_a(val_hap1) == val_member.index_a(val_hap2) and val_member.count_a(val_hap1) == 1):
                    pairs_consistent[idx_allele] = True
                    continue

        if all(pairs_consistent):
            return True
        return False


def is_equal(al_1, al_2):
    """
    check if two strings represent same allele
    """
    if al_1 == al_2 or al_1.startswith(al_2) or al_2.startswith(al_1) or al_1 == '' or al_2 == '':
        return True
    return False


def errors_summary_and_writing_to_file(output_path, errors_in_families, families_num):
    """
    write details of all the errors in families, and summary, to an error file
    :param output_path: output path
    :param errors_in_families:  dict with errors
    :param families_num: number of families in the data
    """
    # cur_path = os.path.abspath("GR_code/GG_GRAMM/data")
    cur_path = os.path.dirname(os.path.realpath(__file__)) + '/data'
    with open(os.path.join(cur_path, 'errors_codes_meaning.json')) as errors_json:
        errors_meaning = json.load(errors_json)

    # ordered keys of families with errors
    sorted_keys = sorted(errors_in_families.keys(), key=lambda my_key: int(my_key))

    with open(''.join([output_path, '/errors.txt']), 'w') as out_file:
        for idx_fam in sorted_keys:
            problematic_member, error_code = errors_in_families[idx_fam]
            error_massage = errors_meaning[error_code]
            if problematic_member == 'All':
                out_file.write(''.join(['Family "', str(idx_fam), '": ', error_massage, '\n']))
            else:
                out_file.write(''.join(['Child "', str(problematic_member), '" in family "', str(idx_fam), '": ',
                                        error_massage, '\n']))

        out_file.write('\n-------------------- Errors Summary (Only rejected families) --------------------\n')
        all_errors_codes = [value[1] for value in errors_in_families.values() if value[0] == 'All']
        for key_error, error_value in errors_meaning.items():
            sum_error = all_errors_codes.count(key_error)
            if sum_error > 0:
                out_file.write(''.join(['Error type: "', error_value, '"\n\tInvalid families: ', str(sum_error),
                                        ' families. ', '{:.1f}'.format(sum_error / families_num * 100),
                                        '% from the data.\n']))
        out_file.write('\nTotal errors:\n' + str(len(all_errors_codes)) + ' from ' + str(families_num) + ' families.\n')
