import copy
from .general_aux_funs import empty_hap


def get_inherited_haps(hapF, hapM, associating):
    """
    interpret the list 'associating' to pair of the haplotypes the _child_ inherited, and pair of not inherited.
    for example, if associating = [2, 1], so inherited haps are hapF.hap2, hapM.hap1. (and not inherited are hapF.hap1, hapM.hap2)
    if the associating is unknown, the value in 'associating' is 'None', so if associating = [None, 2] so inherited are
    None and hapM.hap2
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param associating: associating list
    :return: hapF_inherited, hapF_no_inherited, hapM_inherited, hapM_no_inherited
    """
    inheritance = [None, None, None, None]  # order:hapF_inherited, hapF_no_inherited, hapM_inherited, hapM_no_inherited
    for idx_par, hap, associate in zip([0, 2], [hapF, hapM], associating):
        if associate == 1:  # hap1 is inherited
            inheritance[idx_par] = hap.hap1
            inheritance[idx_par + 1] = hap.hap2
        elif associate == 2:  # hap2 is inherited
            inheritance[idx_par] = hap.hap2
            inheritance[idx_par + 1] = hap.hap1
    return inheritance[0], inheritance[1], inheritance[2], inheritance[3]


def get_high_res_allele(str1, str2):
    """
    return allele with high res. if str1 = '01', str2='01:02', return str2
    :param str1: string of first allele
    :param str2: atring of second allele
    :return: allele with high resolution
    """
    if len(str1) >= len(str2):
        return str1
    return str2


def compare_with_connected_hap_and_remove_if_needed(p_als, no_inher_hap, inher_hap_second_par, no_inher_hap_second_par,
                                                    idx_common, allele_name):
    """
    while we insert data to haplotypes, we often insert 2 values to 2 locations, because we can't determine 1
    (for example: child inherited hap1 of F, hap1 of M, and in allele 'B' they have no data, and child has [01, 02],
    so we insert they both ([01, 02]) to the 2 haps)
    then, if we determine in one haplotype some value (01 in F, for example),
    we want that in the second haplotype will be determine the second (02 in M).

    so, in this function, we go over the other haplotypes (in a specific allele) and check if one of them has two
    identical values such the current haplotype. if so, we remove for it the value that we determine in the current hap.

    but - before we do it, we have to check whether there is more than one hap that has the same values
    (for example - either F and M has [01, 02] in 'B' in hap1, but F has [01, 02] in 'B' in hap2, too. so if we
    determine '01' in hap1 of F, we can't know if determine '02' in hap1 of M or in hap2 of F)
    if so, we don't know which of them were inserted together. in this case we don't remove.

    :param p_als: the allele from current parent, that _child_ inherited
    :param no_inher_hap: the haplotype from the current parent, that the _child_ didn't inherit
    :param inher_hap_second_par: the haplotype that _child_ inherited from second parent
    :param no_inher_hap_second_par: the haplotype that _child_ didn't inherit from second parent
    :param idx_common: idx of the common value between parent and _child_ (that we want determine in p_als and remove
    from the connected haplotype)
    :param allele_name: allele name
    """
    no_inherited_hap_in_par = no_inher_hap[allele_name] if no_inher_hap else None
    inherited_second_par = inher_hap_second_par[allele_name] if inher_hap_second_par else None
    no_inherited_second_par = no_inher_hap_second_par[allele_name] if no_inher_hap_second_par else None

    if no_inherited_hap_in_par != inherited_second_par and no_inherited_hap_in_par != no_inherited_second_par and \
            (inherited_second_par != no_inherited_second_par or not inherited_second_par):  # the last is for cases that they both are None

        # first condition (if no_inherited.., elif inherited_second..) for check it's not None
        # second (p_als == ..) for check that they have same values (were inserted together)
        if no_inherited_hap_in_par and p_als == no_inherited_hap_in_par:
            no_inherited_hap_in_par.remove_a(p_als[idx_common])
        elif inherited_second_par and p_als == inherited_second_par:
            inherited_second_par.remove_a(p_als[idx_common])
        elif no_inherited_second_par and p_als == no_inherited_second_par:
            no_inherited_second_par.remove_a(p_als[idx_common])


def M_hap_is_fuller(hapF_inherited, hapM_inherited, alleles_names):
    """
    check if M haplotype is fuller (= more alleles with values) than F hap.
    """
    if not hapF_inherited or not hapM_inherited:  # one of them is None, so the order is not matter
        return False

    if empty_hap(hapF_inherited):
        return True

    count_fuller = 0
    for allele_name in alleles_names:
        Als_F = hapF_inherited[allele_name]
        Als_M = hapM_inherited[allele_name]
        if (not any(Als_F)) and any(Als_M):  # using 'any' because could be [''] that we want will be consider as empty
            count_fuller += 1
        elif any(Als_F) and (not any(Als_M)):
            count_fuller -= 1
    return True if count_fuller > 0 else False


def add_child_data(hapF, hapM, child, idx_child, children_num, associating, alleles_names, aux_tools, idx_fam, errors_in_families):
    """
    after we associated a _child_ with the haplotypes he inherited from parents, we go over on the alleles_names of _child_ and
    parents in these haplotypes and try to add data to the parents, based on the comparison between them.
    we handle the options :
        empty _child_ (no data in allele): do nothing
        empty parent (no data in allele): add _child_ data to parent
        _child_ and parent both have 1 value in allele: if contradictory - error, else: insert the value with high resolution
        _child_ has 2 parent has 1: if contradictory - error, else: insert the value with high res, and remove from _child_
        the common allele
        _child_ has 1 parent has 2: if contradictory - error, else: insert the value with high res, and remove the non-common
        from parent (+ maybe remove from another haplotype, explained in 'compare_with..._remove_if_needed')
        _child_ has 2 parent has 2: like above, with remove the common allele from _child_

    about errors: if only one _child_ is problematic - continue, if more than one - stop the running.
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param child: _child_ dict
    :param idx_child: _child_ idx (1, 2..)
    :param children_num: number of children in family
    :param associating: associating list
    :param alleles_names: alleles_names names
    :param aux_tools: dict with auxiliary tools
    :param idx_fam: family count (1, 2..)
    :param errors_in_families: dict contains errors in families
    :return: if adding was succeed or not
    """
    _child_ = copy.deepcopy(child)  # copy because we change '_child_'
    hapF_inherited, hapF_no_inherited, hapM_inherited, hapM_no_inherited = get_inherited_haps(hapF, hapM, associating)

    # if 'hapM_inherited' contains more allele than F, we want to begin with M
    # (it's better to begin with fuller haplotype, because in the comparisons, values could be removed from child data,
    # and when we arrive to the more empty hap, the insertion is more efficient (insert 1 value instead of 2))
    if M_hap_is_fuller(hapF_inherited, hapM_inherited, alleles_names):
        inherited_haps, no_inherited_haps = [hapM_inherited, hapF_inherited], [hapM_no_inherited, hapF_no_inherited]
    else:
        inherited_haps, no_inherited_haps = [hapF_inherited, hapM_inherited], [hapF_no_inherited, hapM_no_inherited]

    # inherited_haps, no_inherited_haps = [hapF_inherited, hapM_inherited], [hapF_no_inherited, hapM_no_inherited]
    continue_running = True

    def child_error():
        # if more than one problematic _child_ exist (the value in the dict != None),
        # or only one _child_ exist in family (so we can't allow he will be problematic)
        if aux_tools['problematic_child'] or children_num == 1:
            return False  # stop the running
        else:
            aux_tools['problematic_child'] = idx_child  # this is the first problematic _child_, so we allow it
            return True  # continue the running

    # we run also on the haplotypes from the second parent for calling to function 'compare_with..remove_if_needed'
    for inher_hap, no_inher_hap, inher_hap_second_par, no_inher_hap_second_par in zip(inherited_haps, no_inherited_haps,
                                                                                      inherited_haps[::-1],
                                                                                      no_inherited_haps[::-1]):

        if not inher_hap:  # hap = None, we have no data about the inheritance of this parent
            continue
        for allele_name in alleles_names:
            c_als = _child_[allele_name]
            p_als = inher_hap[allele_name]

            empty_child = True if c_als.empty_Als() else False  # _child_: ['', '']
            empty_parent = True if p_als.empty_Als() else False  # parent: ['', '']
            child_has_1_parent_has_1 = True if len(c_als) == len(p_als) == 1 else False  # c:[01] p:[01:05]
            child_has_2_parent_has_1 = True if (len(c_als) == 2 and len(p_als) == 1) else False   # c:[01:05, 02] p:[01]
            child_has_1_parent_has_2 = True if (len(c_als) == 1 and len(p_als) == 2) else False  # c:[01] p:[01:05, 02]
            child_has_2_parent_has_2 = True if len(c_als) == len(p_als) == 2 else False  # c:[01:05, 02] p:[01, 03]

            if empty_child:
                continue  # to next allele

            elif empty_parent:  # tested
                # if _child_ homozygous, remove one value and insert only one to the parent
                if len(c_als) == 2 and c_als[0] == c_als[1]:
                    c_als.remove_a(c_als[1])
                p_als.clear()  # p_als: ['', ''] - > []. for not contain empty value when we add new values
                p_als.extend(c_als)  # insert _child_ data to parent

            elif child_has_1_parent_has_1:
                if c_als != p_als:  # the values are contradictory, so it's an error
                    continue_running = child_error()
                else:
                    p_als[0] = get_high_res_allele(p_als[0], c_als[0])  # insert to parent the value with high res

            elif child_has_2_parent_has_1:
                if p_als[0] not in c_als:  # the values are contradictory, so it's an error
                    continue_running = child_error()
                else:
                    idx_common = c_als.index_a(p_als[0])
                    p_als[0] = get_high_res_allele(p_als[0], c_als[idx_common])
                    c_als.remove_a(c_als[idx_common]) # remove common allele from _child_- because it's inserted to parent

            elif child_has_1_parent_has_2:  # tested
                if c_als[0] not in p_als:
                    continue_running = child_error()  # the values are contradictory, so it's an error
                else:
                    idx_common = p_als.index_a(c_als[0])
                    # insert the value with the high res
                    p_als[idx_common] = get_high_res_allele(p_als[idx_common], c_als[0])
                    compare_with_connected_hap_and_remove_if_needed(p_als, no_inher_hap, inher_hap_second_par,
                                                                    no_inher_hap_second_par, idx_common, allele_name)
                    p_als.remove_a(p_als[1 - idx_common])  # remove the non-common allele between parent and _child_

            elif child_has_2_parent_has_2:  # tested
                par_inter, idx_par_inter = c_als.intersection(p_als)
                if par_inter == 1:  # one intersection between values of _child_ and parent
                    idx_common_child = c_als.index_a(p_als[idx_par_inter])
                    p_als[idx_par_inter] = get_high_res_allele(p_als[idx_par_inter], c_als[idx_common_child])

                    compare_with_connected_hap_and_remove_if_needed(p_als, no_inher_hap, inher_hap_second_par,
                                                                    no_inher_hap_second_par, idx_par_inter, allele_name)
                    p_als.remove_a(p_als[1 - idx_par_inter])  # remove the non-common allele between parent and _child_
                    c_als.remove_a(c_als[idx_common_child])  # remove the common allele between parent and _child_

                elif par_inter == 0:
                    continue_running = child_error()

            if not continue_running:  # more than 1 problematic _child_ exist
                break

        if not continue_running:  # more than 1 problematic _child_ exist
            break

    success = True if continue_running else False  # unnecessary variable, but I added for the readability

    if not success:
        errors_in_families[idx_fam] = ['All', '7']

    return success







