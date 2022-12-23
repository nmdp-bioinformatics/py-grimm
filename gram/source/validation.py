from .als_update import Als


def is_valid_family(family, alleles_names, par_num, count_fam, aux_tools, errors_in_families):
    """
    call to 'data_validation()' that check the family validation.
    if the family is valid - return True
    else: check if there are only one problematic _child_.
        if it's only one _child_: remove him and continue the analyzing. (return True)
        but if the problem is with a parent or some children - reject the family (return False)
    :param family: family dict
    :param alleles_names: alleles_names names
    :param par_num: num of parents
    :param count_fam: family index
    :param aux_tools: dict that we can sign if there is a problematic _child_
    :param errors_in_families: dict that includes the details about errors
    :return: True if valid, otherwise - False
    """
    is_valid, error_code, invalid_member = data_validation(family, alleles_names, par_num)
    if not is_valid:
        # if the error relates to a parent, or more than one _child_, reject the family
        if invalid_member == 'All':
            errors_in_families[count_fam] = ['All', error_code]
            return False
        # if there is only one problematic _child_, we mark him and remove temporarily from family.
        # (in the output to user, we will sign this _child_ as problematic)
        else:
            errors_in_families[count_fam] = [invalid_member, error_code]
            del family[invalid_member]  # it is removed also from families_dict (because of 'del' attribute)
            aux_tools['problematic_child'] = invalid_member
    return True


def data_validation(fam_d, alleles_names, par_num):
    """
    check validation of the user data.
    there are some kinds of errors, signed with errors codes. the meaning of them is in "errors_codes_meaning.json"
    :param fam_d: family dict
    :param alleles_names: alleles_names names list: [A, B, C, DRB1, DQB1]
    :param par_num: number of parents
    :return: 3 outputs:
        first:
            True if id, False otherwise
        second:
            the error code ('None' if valid)
        third:
            the family member that cause to the error. ('None' if valid)
            if there is only one problematic _child_, we want to analyze
            the family without him.
            in any other case (errors in parents, or errors that relate for some members), we reject the family,
            and return 'All' as a flag
    """
    invalid_cases = []  # in this list we save all the invalid cases in this family

    # validation tests. the order is important
    check_invalid_character(fam_d, alleles_names, invalid_cases)
    check_missing_data(fam_d, par_num, invalid_cases)
    check_too_much_alleles(fam_d, alleles_names, invalid_cases)

    if par_num == 2:
        check_allele_in_child_that_does_not_exist_in_parents(fam_d, alleles_names, par_num, invalid_cases)

    if par_num == 0:
        check_if_there_is_allele_with_4_diff_values(fam_d, alleles_names, invalid_cases)

    if len(invalid_cases) == 0:
        return True, None, None

    invalid_members = [case[1] for case in invalid_cases]  # get the invalid family members (e.g: ['F', '2'])

    # in case that there is an error only in one _child_,
    # we want to sign and remove him, and analyze the family without him
    if set(invalid_members) == 1 and invalid_members[0] not in ['F', 'M', 'All']:
        # return that it's invalid, the error code, and the _child_ index
        return False, invalid_cases[0][0], invalid_members[0]
    # in any other case, we want to reject the family.
    # if there are some reasons for the rejection, we return the first
    return False, invalid_cases[0][0], 'All'


def check_invalid_character(fam_d, alleles_names, invalid_cases):
    """
    check if there is an invalid character in data
    """
    for al_name in alleles_names:  # [A, B..]
        for fam_member in fam_d:  # [F, M, 1..]
            for single_al in fam_d[fam_member][al_name]:  # [02:01, 30:04]
                if single_al == "" or single_al == " ":
                    continue
                elif ":" in single_al:
                    parts = single_al.split(":")
                    for part in parts:
                        if not (part.isnumeric() or part.isupper() or part == ''):
                            invalid_cases.append(('1', fam_member))  # invalid character in data
                else:
                    if not (single_al.isnumeric() or single_al.isupper()):
                        invalid_cases.append(('1', fam_member))  # invalid character in data


def check_missing_data(fam_d, par_num, invalid_cases):
    """
    check for missing data
    """
    if len(fam_d) < 2:
        invalid_cases.append(('2', 'All'))  # less than two people
    if len(fam_d) == 2 and par_num == 2:
        invalid_cases.append(('3', 'All'))  # no children


def check_too_much_alleles(fam_d, alleles_names, invalid_cases):
    """
    check if there are too much alleles_names in the family (more than 4 in an allele)
    """
    for al_name in alleles_names:
        lst = Als()
        for fam_member in fam_d:  # [F, M, 1 ...]
            if any(fam_d[fam_member][al_name]):  # not empty
                # lst = fam_d[fam_member][al_name].merge(lst)
                lst = lst.merge(fam_d[fam_member][al_name])
        if len(lst) > 4:
            invalid_cases.append(('4', 'All'))  # Too many alleles_names


def check_allele_in_child_that_does_not_exist_in_parents(fam_d, alleles_names, par_num, invalid_cases):
    """
    check if there is an allele in a _child_ that does not exist in the parents
    """
    for al_name in alleles_names:
        fm_als = fam_d['F'][al_name] + fam_d['M'][al_name]
        for fam_member in fam_d:
            if fam_member != 'F' and fam_member != 'M' and len(fm_als) == 4 and all(fm_als):
                in_fm = fam_d[fam_member][al_name].sub_lst(fm_als)
                if not in_fm:
                    invalid_cases.append(('5', fam_member))  # allele in a _child_ that does not exist in the parents


def check_if_there_is_allele_with_4_diff_values(fam_d, alleles_names, invalid_cases):
    """
    if there are no parents, and there isn't an allele with 4 different values in the children,
    the algorithm could not be executed
    (because when no parents, we rely on the assumption that there are 4 different in one allele, at least)
    it's not an "invalid" input, but we reject the family because we can not analyze it
    """
    four_different_values = False  # flag to know if there is allele with 4 different values

    for al_name in alleles_names:
        if four_different_values:  # if there is allele with 4 diff, do not need to check more
            return

        alleles_values = Als()
        for child in fam_d:  # no parents, according the condition in the call to this function
            child_alleles = fam_d[child][al_name]
            if any(child_alleles):  # merge only if there is data in 'child_alleles'
                # merge the values of the current alleles_names _child_ to the values of the other children in this allele
                # e.g. : alleles_values: [02, 03]. child_alleles: [02:01, 04]. so the merging: [02, 03, 04]
                # note: the merging may save the low-res values (02 instead of 02:01),
                # but it's not matter, because just need the values amount
                alleles_values = alleles_values.merge(child_alleles)
        # after go over on all the children values in the current allele, check if there are 4 values
        if len(alleles_values) == 4:
            four_different_values = True

    if not four_different_values:
        invalid_cases.append(('6', 'All'))  # no parents, and no alleles_names with 4 diff values (algorithm can not executed)

