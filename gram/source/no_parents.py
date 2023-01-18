from .als_update import Als

"""
Calls to function order (2 clusters of functions):

cluster 1:
-> insert_data_to_single_allele_in_par_haps (from runfile)
    -> create_merged_allele_values
    -> check_num_child_with_diff_values_in_specific_allele
    -> insert_data_of_2_children_no_homozygous
        -> get_4_alleles_from_2_children
    -> insert_data_of_2_children_one_homozygous
        -> get_4_alleles_from_2_children
    -> insert_data_more_than_2_children
        -> divide_alleles_to_2_groups
            -> check_if_exist_homoz
            -> divide_2_alleles_to_non_empty_groups

cluster 2:                   
-> associate_children_to_par_haps_while_parents_dont_exist (from runfile)
    ->rest_of_haps_is_empty
"""


def create_merged_allele_values(children, al_name):
    """
    create a list (Als) with all the different value in children data about specific allele
    for example: if _child_ 1 has A : [02:01, 03:04], _child_ 2: A: [05, 03], _child_ 3: A: [02, 07],
    so lst will be [02:01, 03:04, 05, 07]
    :param children: children dict
    :param al_name: name of specific allele
    """
    lst = Als()
    for child in children:
        lst = lst.merge(children[child][al_name])
    return lst


def check_num_child_with_diff_values_in_specific_allele(children, al_name):
    """
    for a specific allele (that we know there are 4 diff values about it in children data),
    we want to now how many children have different values in this allele
    for example, _child_ 1: A[01, 02], _child_ 2: A[03, 04], _child_ 3: A[01, 02]
    so the num is 2 children (cause 1 and 3 are identical in the values in this allele)
    :param children: children dict
    :param al_name: name of specific allele
    :return: children dict, that contains data about one allele only
    """
    exist = []
    dict_children_one_allele = {}
    for child in children:
        # check that the values are not empty (["", ""]), and not in 'exist'
        if any(children[child][al_name]) and children[child][al_name] not in exist:
            dict_children_one_allele[child] = children[child][al_name]
        exist.append(children[child][al_name])
    return dict_children_one_allele


def get_4_alleles_from_2_children(dict_children_one_allele):
    """
    for example: dict = {'1': [01, 02], '2':[03, 04]}, so return 01, 02, 03, 04
    :param dict_children_one_allele: children dict, that contains data about one allele only
    """
    alleles = []
    for alleles_child in dict_children_one_allele.values():
        alleles.extend(alleles_child)
    return alleles[0], alleles[1], alleles[2], alleles[3]


def insert_data_of_2_children_no_homozygous(dict_children_one_allele, hapF, hapM, al_name):
    """
    when we have 2 children with 4 different values in an allele, we insert the data to parents haplotypes such that:
    _child_ 1 : [X1, X2], _child_ 2: [Y1, Y2] --> hapF.hap1 = [X1], hapM.hap1 = [X2]; hapF.hap2 = hapM.hap2 = [Y1, Y2]
    :param dict_children_one_allele: children dict, that contains data about one allele only
    :param hapF: haplotype father
    :param hapM: haplotype mother
    :param al_name: allele name
    """

    child1_al1, child1_al2, child2_al1, child2_al2 = get_4_alleles_from_2_children(dict_children_one_allele)

    hapF.hap1[al_name].append(child1_al1)
    hapM.hap1[al_name].append(child1_al2)

    hapF.hap2[al_name].extend([child2_al1, child2_al2])
    hapM.hap2[al_name].extend([child2_al1, child2_al2])


def insert_data_of_2_children_one_homozygous(dict_children_one_allele, hapF, hapM, al_name):
    """
    like the function above, just in case that one of the children is homozygous. so the insertion is deterministic.
    for example: _child_ 1: [01, 01], _child_ 2: [02, 03] so the parents are 01+02, 01+03
    """
    child1_al1, child1_al2, child2_al1, child2_al2 = get_4_alleles_from_2_children(dict_children_one_allele)

    for parent_hap, child_allele in zip([hapF.hap1, hapF.hap2, hapM.hap1, hapM.hap2],
                                        [child1_al1, child2_al1, child1_al2, child2_al2]):
        parent_hap[al_name].append(child_allele)


def check_if_exist_homoz(dict_children_one_allele):
    for alleles in dict_children_one_allele.values():
        if alleles[0] == alleles[1]:
            return True, alleles[0]
    return False, None


def divide_2_alleles_to_non_empty_groups(gr1, gr2, alleles, try_again):
    """
    check the intersection between g1, g2 (that have at least one allele in each of them) and the alleles_names.
    if no intersection (g1:[01], g2:[02], alleles_names:[03, 04]: we cannot add data to groups, so we save the alleles_names in the
    list 'try_again', and we call to this function again with those alleles_names (in hope that in the second time, after
    adding more information, will be intersection.)
    if 1 intersection (g1:[01], g2[02], alleles_names:[01, 03]: so we add the allele who not the intersection ('03'), to the
    group that do not have the intersection (g2)
    if 2 intersections: it means that there is not new information to gr1 and gr2
    :param gr1: group 1
    :param gr2: group 2
    :param alleles: alleles_names pair
    :param try_again: list that let us know if we need to execute this function again
    :return:
    """
    gr1_gr2 = gr1 + gr2
    intersection_count, idx_intersection = gr1_gr2.intersection(alleles)

    if intersection_count == 0:
        try_again.append(alleles)
    elif intersection_count == 1:
        common_allele = alleles[idx_intersection]
        not_common_allele = alleles[1 - idx_intersection]
        if common_allele in gr1 and common_allele not in gr2:
            gr2.append(not_common_allele)
        else:
            gr1.append(not_common_allele)


def divide_alleles_to_2_groups(dict_children_one_allele):
    """
    dividing alleles_names of 3 children or more to 2 groups (one for each parent)
    for example, c1:[01, 02], c2:[02, 03], c3:[01, 04], so -> par1: 01~03, par2: 02~04
    there are 2 cases:
        1. easy case: there is homozygous _child_, so divide his alleles_names to the 2 groups, and then go over the other
        children alleles_names, and insert in some order (no matter how) to the groups, until each group is of size 2
        for example: c1:[01, 01], c2:[01, 02], c3:[02, 03]
        -->(iter1) par1: 01, par2: 01  -->(iter2) par1: 01~02, par2: 01  -->(iter3) par1: 01~02, par2: 01~03
        2. difficult case: no homozygous _child_, so add the alleles_names of the first _child_, and then, for the others, call
        to 'divide_2_alleles_to_non_empty_groups'
    :param dict_children_one_allele: children dict, that contains data about one allele only
    :return: 2 groups
    """
    gr1, gr2 = Als(), Als()

    is_homoz, homoz_allele = check_if_exist_homoz(dict_children_one_allele)
    if is_homoz:  # if a _child_ has [01, 01] so each parent has '01'
        gr1.append(homoz_allele)
        gr2.append(homoz_allele)
        for alleles in dict_children_one_allele.values():
            for al in alleles:
                if al not in gr1 and al not in gr2:
                    if len(gr1) < 2:
                        gr1.append(al)
                    elif len(gr2) < 2:
                        gr2.append(al)
            if len(gr1) == len(gr2) == 2:  # the groups are full
                break

    else:  # no homozygous
        # 'try_again' let us know if one allele did not succeed to be inserted to the groups in first time, so after
        # insertion the other, we try to insert it again. more explanations in documentation of 'divide_2_alleles..'
        try_again = []
        for alleles in dict_children_one_allele.values():
            if len(gr1) == len(gr2) == 0:  # first insertion, the order is not matter
                gr1.append(alleles[0])
                gr2.append(alleles[1])
            else:
                divide_2_alleles_to_non_empty_groups(gr1, gr2, alleles, try_again)
        if len(try_again) > 0:
            divide_2_alleles_to_non_empty_groups(gr1, gr2, try_again[0], [])

    return gr1, gr2


def insert_data_more_than_2_children(dict_children_one_allele, hapF, hapM, al_name):
    """
    when we have 3 or more children, we try to divide their allele values to 2 groups:
    one that were inherited from parent 1, and second from parent 2
    :param dict_children_one_allele: children dict, that contains data about one allele only
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param al_name: name of the allele in 'dict_children_one_allele'
    :return: True if we success to divide, False otherwise
    """
    group1, group2 = divide_alleles_to_2_groups(dict_children_one_allele)

    if len(group1) < 2 or len(group2) < 2:
        return False

    # insert data to parents haplotypes
    for parent_hap, allele_from_group in zip([hapF.hap1, hapF.hap2, hapM.hap1, hapM.hap2],
                                             [group1[0], group1[1], group2[0], group2[1]]):
        parent_hap[al_name].append(allele_from_group)
    return True


def insert_data_to_single_allele_in_par_haps(family, alleles_names, hapF, hapM, count_fam, errors_in_families):
    """
    in case of no parents, we try to found an allele that children have 4 different values.
    (for example, in 'B', children data includes: 01, 02, 03, 04)
    thus, we can (in some terms) insert this data to parents haplotypes
    :param family: family dict
    :param alleles_names: alleles_names names
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param count_fam: family index
    :param errors_in_families: errors dict, for report if we dont find 4 different values
    """

    success = False  # flag to know if we succeeded
    allele_inserted_to_parents = None  # the name of the allele that we insert the data (if we succeeded)
    # is_deter = True  # flag that =False if there are 2 children with 4 diff values,  #TODO: remove 'is_deter'
    # # so the insertion is not deterministic (par1: [al1_child1, al1/al2_child2], par2: [al2_child1, al1/al2_child2])
    # # it's important for a function later, that associate children with parents haplotypes

    children_keys = [k for k in family.keys() if k not in ['F', 'M']]
    children = {k: family[k] for k in children_keys}  # get dub-dict of 'family', only with children

    for al_name in alleles_names:
        merged_allele_values = create_merged_allele_values(children, al_name)
        if len(merged_allele_values) == 4:

            # 'dict_children_one_allele' includes children data about single allele, with no repetitions
            # and it does not includes children with empty data in this allele (["", ""])
            dict_children_one_allele = check_num_child_with_diff_values_in_specific_allele(children, al_name)
            num_children = len(dict_children_one_allele)

            # when num_children = 2, there are 2 options: two children heterozygous,
            # or an homozygous _child_ (could be also the 2 children).
            # in the first case, the insertion is not deterministic (explained in the function 'insert..no_homozygous')
            # in the second, the insertion is deterministic (explained in the function 'insert..one_homozygous')
            if num_children == 2:
                if len(set(merged_allele_values)) == 4:  # completely different
                    insert_data_of_2_children_no_homozygous(dict_children_one_allele, hapF, hapM, al_name)
                elif len(set(merged_allele_values)) < 4:  # there is homozygous _child_
                    insert_data_of_2_children_one_homozygous(dict_children_one_allele, hapF, hapM, al_name)
                success = True
                allele_inserted_to_parents = al_name
                break  # we do not need to check more alleles_names, if we succeeded.

            else:  # num_children > 2. 1 is not an option, because he can't have 4 diff values in an allele.
                success = insert_data_more_than_2_children(dict_children_one_allele, hapF, hapM, al_name)
                if success:
                    allele_inserted_to_parents = al_name
                    break  # we do not need to check more alleles_names, if we succeeded.

    if not success:
        errors_in_families[count_fam] = ['All', '6']

    return success, allele_inserted_to_parents


def rest_of_haps_is_empty(haps, allele_inserted):
    """
    go over parent haplotype and check if they are empty, except the allele of the first insertion
    (with children data)
    :param haps: dual haplotype of a parent
    :param allele_inserted: an allele with data, that we do not want to check
    :return: True if the rest is empty, False otherwise
    """
    for hap in [haps.hap1, haps.hap2]:
        for allele_name, allele_data in hap.items():
            if allele_name != allele_inserted and allele_data:
                return False
    return True


def associate_children_to_par_haps_while_parents_dont_exist(hapF, hapM, child, allele_inserted):
    """
    associate _child_ to the parents haplotypes he inherited. this function is for scenario of no parents in data.
    in this scenario, there are 1 allele with data which were inserted before ('allele_inserted'),
    and we try to associate according to him.
    by list 'associating' we mark the haplotypes the _child_ inherited (in first index: father hap, in second: mother)

    we handle with some cases:
    1. _child_ with no data about this allele (["", ""]) - cannot be associated
    2. each value in _child_ matches to a value in parents haps (f:01~02 m:03~04, c:[01,04], so associating: [1, 2])
    3. each value in _child_ matches to 2 values in parents haps. it happened in specific scenario:
       in previous stage were 2 children with 4 different values, so we inserted them non-deterministically
       (if c1:[01, 02], c2:[03, 04], so par1: 01~03/04, par2: 02~03/04, so if we associate c2, it's [2, 2])
    4. one value in _child_ has 1 match, and the second: 2 matches. it could happened in one of the following options:
        4a. one parent is homozygous (f:01~01, m:02~03, c:[01, 02]) so we not sure which hap (in the homoz parent) we
            need associate with the _child_ (in f: 1 or 2?).
            in this case, we check if it's the first _child_ in the external loop, (by check if the rest of the haplotypes
            of this parent are empty). if so, we associate randomly (in index 1), because it doesnt matter.
            but otherwise, we do not associate.
        4b. there is common value between the parents, and _child_ has this value (f:01~02, m:02~03, c:[01, 02]),
            so first we associate with the parent with 1 match (m, in the example). then, we take the non-common
            allele (01) and check its index in the parent with the 2 matches.
    * other cases look impossible, but maybe I missed something. anyway, if we don't meet one of the conditions
      above, there is no associating.
      (about the cases: 2 parents homozygous (so 2 matches twice); or _child_ + 1 par are homoz' (so 1 match and 2 matches))
      they shouldn't appear, because according the previous stage ('insert_data_to_single...') ,their insertion
      couldn't happened.
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param child: _child_ data
    :param allele_inserted: the allele with the data in parents haplotypes
    :return: True if the associating succeed, False otherwise ; and 'associating' list
    """
    associating = [None, None]  # first index: father, second: mother

    child_alleles = child[allele_inserted]  # _child_ values of 'allele_inserted'
    f_haps = hapF.hap1[allele_inserted] + hapF.hap2[allele_inserted]
    m_haps = hapM.hap1[allele_inserted] + hapM.hap2[allele_inserted]

    # check intersections between _child_ alleles_names and parents haps values
    f_inter, idx_f_inter = child_alleles.intersection(f_haps)
    m_inter, idx_m_inter = child_alleles.intersection(m_haps)

    # case 1
    if not any(child_alleles):  # _child_ has no data in this allele
        pass  # this condition is just for not check the other conditions

    # case 2
    elif f_inter == m_inter == 1:
        associating[0] = idx_f_inter + 1
        associating[1] = idx_m_inter + 1

    # case 3
    elif f_inter == m_inter == 2:
        associating[0] = 2
        associating[1] = 2

    # case 4
    elif (f_inter == 1 and m_inter == 2) or (f_inter == 2 and m_inter == 1):
        # 4a
        if f_haps[0] == f_haps[1] or m_haps[0] == m_haps[1]:
            haps_par_homoz = hapF if f_inter > m_inter else hapM
            if rest_of_haps_is_empty(haps_par_homoz, allele_inserted):
                if f_inter > m_inter:  # f_inter == 2
                    associating[0] = 1
                    associating[1] = idx_m_inter + 1
                else:  # m_inter == 2
                    associating[0] = idx_f_inter + 1
                    associating[1] = 1
        # 4b
        else:
            if f_inter < m_inter:  # f_inter == 1
                associating[0] = idx_f_inter + 1
                associating[1] = 1 if m_haps[0] not in f_haps else 2
            else:  # m_inter == 1
                associating[0] = 1 if f_haps[0] not in m_haps else 2
                associating[1] = idx_m_inter + 1

    success = True if all(associating) else False

    return success, associating

