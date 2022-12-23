import copy
import itertools
import pyard
from .general_aux_funs import empty_hap
from .als_update import Als


def duplicate_hap_if_one_empty(hapF, hapM, aux_tools, count_fam):
    """
    GRIMM cannot handle with case of one haplotype with data and one empty.
    but this case is possible, when we have 1 parent and 1 _child_, for example. the haplotype that the _child_ didn't
    inherit from the second parent, has no data.
    therefore, the solution is to duplicate the haplotype with the data, and after GRIMM return to original condition
    (clear one haplotype). it's is handled by the part of 'GG_Post'
    :param hapF: father haplotypes
    :param hapM: mother haplotypes
    :param aux_tools: dict with auxiliary tools
    :param count_fam: idx of the family
    """
    h1F, h2F, h1M, h2M = hapF.hap1, hapF.hap2, hapM.hap1, hapM.hap2
    for idx_par, (hap1, hap2), dual_hap in zip(['F', 'M'], [(h1F, h2F), (h1M, h2M)], [hapF, hapM]):
        if empty_hap(hap1):
            dual_hap.hap1 = copy.deepcopy(hap2)
            aux_tools['parent_has_empty_hap'][count_fam] = idx_par  # sign it for treatment later
        elif empty_hap(hap2):
            dual_hap.hap2 = copy.deepcopy(hap1)
            aux_tools['parent_has_empty_hap'][count_fam] = idx_par  # sign it for treatment later


def remove_data_if_just_one_allele_full(h1, h2):
    """
    GRIMM cannot handle with case of allele that have data in one haplotype but no data in the second haplotype
    for example: A*02:01+A*ZZZZ (ZZZZ means empty).
    so in alleles_names like this, we delete the data from the allele in the full haplotype
    - pay attention: this function is called after the function 'duplicate_hap_if_one_empty', because otherwise,
      we might delete akk the data from the haplotypes
    :param h1: first haplotype
    :param h2: second haplotype
    """
    for (key1, value1), (key2, value2) in zip(h1.items(), h2.items()):
        if value1.empty_Als() and not value2.empty_Als():
            h2[key2] = Als()
        elif value2.empty_Als() and not value1.empty_Als():
            h1[key1] = Als()


def create_binary_to_GRIMM(h1, h2, alleles_names, par_num):
    """
    we send to GRIMM a binary list (consist of 0/1) for each parent. that because of this case:
    when in some allele, there are 2 identical values in the haplotypes of a parent, e.g:
    hap1: A: [01], B: [03, 04] ...
    hap2: B: [02], B:[03, 04] ...
    so if we separate to 4 options (hap1:01^03/01^04, hap2:02^03/02^04), we get that the first option of hap1 and the
    second of hap2 are same (complement to the same case, because if 01 go with 03, so 02 go with 04), and same for
    second of hap1 and first of hap2. so actually we have 2 options, but we don't have a way to show that by the gl-string.
    so in these cases, we sign '1' in the index of the allele with 2 (identical) values, and later we will remove one
    value from each hap (04 from hap1 and 03 from hap2), and when GRIMM will get the sign '1' it will know to open
    the two options (01^03, 02^04 or 01^04, 02^03)
    - there is am extreme case: if there are no parents, and all children have 2 same values in some allele, for example
    03,05 in 'C', so in all parents haps we inserted '03/05' but we cant know if each parent is heterozygous (03+05)
    or homozygous (03+03) so in this case, we leave it with the two options, and GRIMM will open all the options.
    :param h1: first haplotype
    :param h2: second haplotype
    :param alleles_names: alleles_names names
    :param par_num: number of parents
    :return: binary list
    """
    binary = [0] * len(alleles_names)
    for idx_allele, allele_name in enumerate(alleles_names):
        # 2 identical values in alleles_names from the two haplotypes
        # and there are parents in the data (explained in documentation)
        if len(h1[allele_name]) == 2 and h1[allele_name] == h2[allele_name] and par_num != 0:
            binary[idx_allele] = 1
    binary[3], binary[4] = binary[4], binary[3]  # # replace 2 last because in GG_GRIMM, DQ before DR
    return binary


def check_if_one_hap_was_empty(aux_tools, count_fam, idx_par):
    idx_par = 'F' if idx_par == 0 else 'M'
    if count_fam in aux_tools['parent_has_empty_hap']:
        if aux_tools['parent_has_empty_hap'][count_fam] == idx_par:
            return True
    return False


def gl_cartesian_product(h1, h2, _alleles_names_):
    """
    from two haplotypes, create all options to a gl string (cartesian product)
    e.g. (with 3 alleles_names, for simplification) hap1: A:[01, 03], B:[02], C:[04]
                                              hap2: A:[06], B:[07], C:[04, 05]
         --> 'pairs_in_allele_from_2_haps' of allele_name = A: [01+06, 03+06] (same idea for B, C)
         --> all_pairs_all_alleles: [[01+06, 03+06], [02+07], [04+04, 04+05]]
         --> cartesian_options: [[01+06, 02+07, 04+04], [01+06, 02+07, 04+05], [03+06, 02+07, 04+04], [03+06, 02+07, 04+05]]
    :param h1: first haplotype
    :param h2: second haplotype
    :param _alleles_names_: list of alleles_names names, without the allele with no data in two haplotypes
    :return: all options to gl string
    """
    all_pairs_all_alleles = []
    cartesian_options = []
    for allele_name in _alleles_names_:
        pairs_in_allele_from_2_haps = []
        for al1 in h1[allele_name]:
            for al2 in h2[allele_name]:
                pairs_in_allele_from_2_haps.append(''.join([str(al1), '+', str(al2)]))
        all_pairs_all_alleles.append(pairs_in_allele_from_2_haps)
    for option in itertools.product(*all_pairs_all_alleles):
        cartesian_options.append(option)

    return cartesian_options


def gl_cartesian_product_1empty(h1):
    """
    like 'gl_cartesian_product', just in case that was 1 empty haplotype, and we duplicated it.
    so in this case, we handle with duplicate haplotype instead of 2 different.
    e.g.: hap1: A:[01, 03], B:[02], C:[04, 05]
          hap_list = [[01, 03], [02], [04, 05]]
          option (one example) = [01, 02, 04] -> [01+01, 02+02, 04+04]
          cartesian_options = [[01+01, 02+02, 04+04], [01+01, 02+02, 05+05], [03+03, 02+02, 04+04], [03+03, 02+02, 05+05]]
    :param h1: haplotype (that was duplicated)
    :return: all options to gl string
    """
    cartesian_options = []
    hap_list = [list(alleles) for alleles in h1.values()]
    for option in itertools.product(*hap_list):
        option = [''.join([str(op), '+', str(op)]) for op in option]
        cartesian_options.append(option)
    return cartesian_options


def options_to_gl(options, allele_start):
    """
    construct one allele sequence from first 2 digits + many options to second 2 digits
    A*02 + [01, 05] --> A*02:01/A*02:05
    :param options: options to second 2 digits
    :param allele_start: first 2 digits
    :return: full sequence
    """
    if not allele_start.endswith(':'):
        allele_start = ''.join([allele_start, ':'])
    gl = ''.join([allele_start, options[0]])
    if len(options) == 1:
        return gl
    for opt in options[1:]:
        gl += ''.join(['/', allele_start, opt])
    return gl


def single2high(allele, alleles_name, idx_fam, idx_par, open_ambiguity_sim, aux_tools, ard):
    """
    convert allele value to high res, if needed
    :param allele: allele value
    :param alleles_name: allele name
    :param idx_fam: family index
    :param idx_par: parent index (F/M)
    :param open_ambiguity_sim: dict with tho options to ambiguity. used in simulations HARD, HARDER, HARDEST
    :param aux_tools: dict with auxiliary tools
    :param ard: pyard object
    :return: converted allele value
    """
    low2high_dict = aux_tools['low2high']
    amb_dict = aux_tools['amb']
    new_allele = ''
    options = []
    success = True
    idx_par = 'F' if idx_par == 0 else 'M'

    #  ---- allele in high res, e.g: 01:02, so no processing required ----
    if allele.count(':') == 1 and not allele.endswith(':'):
        new_allele = ''.join([alleles_name, '*', allele])
    #  -------------------------- high res --------------------------/

    # ---- allele in low res, e.g. 02. check options in open_ambiguity_sim or pyard or low2high_dict ----
    elif ':' not in allele:
        # one of the simulations, that we removed all the open ambiguities (03:01/03:02.. ->03)
        if open_ambiguity_sim:
            # parent was in simulation file
            if idx_par in open_ambiguity_sim[idx_fam] and allele in open_ambiguity_sim[idx_fam][idx_par][alleles_name]:
                options = open_ambiguity_sim[idx_fam][idx_par][alleles_name][allele]
            # parent wasn't in simulation file but there is data from children genotypes
            # (or exists, but no data in this allele (options=[]))
            elif (idx_par not in open_ambiguity_sim[idx_fam] or not options) and allele in \
                    open_ambiguity_sim[idx_fam]['children'][alleles_name]:
                options = open_ambiguity_sim[idx_fam]['children'][alleles_name][allele]

        if not open_ambiguity_sim or not options:  # not a simulation or a simulation without options to this allele
            try:  # py-ard
                tmp_allele_name = copy.copy(alleles_name)[:2] if alleles_name in ['DRB1', 'DQB1'] else alleles_name
                key_serology = ''.join([tmp_allele_name, copy.copy(allele).lstrip('0')])
                new_allele = ard.redux_gl(key_serology, 'lgx')
            except:  # low2high
                key = ''.join([alleles_name, '*', allele])
                if key in low2high_dict:
                    options = low2high_dict[key]
                else:  # if not in the dict, the allele stays with low res
                    new_allele = key
    # -------------------------- low res --------------------------/

    # ---------- ambiguity, e.g: 02:APC (in the code its 02:, because we removed the letters) ----------
    elif allele.endswith(':'):
        # 'amb_dict' contains the ambiguity of alleles_names in current family that removed before.
        # e.g: in data: A*02:APC , after remove: 02: , in 'amb_dict' (to family 1, for example): {'1': {'A*02': APC}}
        key = ''.join([alleles_name, '*', copy.copy(allele).rstrip(':')])
        ambiguity = amb_dict[idx_fam][key]

        try:  # py-ard
            pyard_key = ''.join([key, ':', ambiguity])
            new_allele = ard.redux_gl(pyard_key, 'lgx')
        except:
            # if ambiguity not in low2high_dict (/pyard?), we open all the options by the low res
            # (e.g: we have A*02:APC but APC not in low2high_dict, so we look for all options of 02--> 02:01/03....)
            if key in low2high_dict:
                options = low2high_dict[key]
            else:
                new_allele = key
    # ------------------------------------ ambiguity ------------------------------------/

    else:
        success = False

    if not success:
        return success, alleles_name

    if new_allele != '':  # we don't need to open options
        return success, new_allele
    else:
        allele_start = ''.join([alleles_name, '*', allele])
        if len(options) >= 1:
            new_allele = options_to_gl(options, allele_start)
        return success, new_allele


def create_gl_string(h1, h2, alleles_names, idx_fam, idx_par, par_num, open_ambiguity_sim, aux_tools, ard):
    """
    from 2 haplotypes (of a parent) we create all the possible gl strings
    :param h1: first haplotype
    :param h2: second haplotype
    :param alleles_names: alleles_names names
    :param idx_fam: family index
    :param idx_par: parent index (0 to father, 1 to mother)
    :param par_num: number of parents
    :param aux_tools: dict with auxiliary tools
    :param open_ambiguity_sim: dict with tho options to ambiguity. used in simulations HARD, HARDER, HARDEST
    :return: if we succeed, and the options to gl strings
    """

    success_gl_string = True
    one_hap_was_empty = check_if_one_hap_was_empty(aux_tools, idx_fam, idx_par)

    # -------------------- open haplotypes to all the options of allele combinations --------------------
    _alleles_names_ = copy.deepcopy(alleles_names)  # copy for remove from it the alleles_names without data (and use later)
    for allele_name in alleles_names:
        if h1[allele_name].empty_Als() and h2[allele_name].empty_Als():  # allele in two haps is empty
            del h1[allele_name]
            del h2[allele_name]
            _alleles_names_.remove(allele_name)
        elif len(h1[allele_name]) == 2 and h1[allele_name] == h2[allele_name] and not one_hap_was_empty and par_num != 0:
            # when 2 same allele in 2 haplotypes (h1: A:[01, 02], h2: A:[01, 02])
            # we determine 1 value in each allele (--> h1: A:[01], h2: A:[02])
            # and in the function 'create_binary_to_GRIMM' sign it by 1, for sign to GRIMM to open this phase.
            # note: we check 'not one_hap_was_empty' because if one was empty, we duplicated and we dont want to remove
            # and check that parents exist in data (otherwise, we dont remove. explained in 'create_binary_to_GRIMM')
            second_allele = h1[allele_name][1]
            h2[allele_name].clear()
            h2[allele_name].append(second_allele)
            h1[allele_name].remove_equal(second_allele)

    if not one_hap_was_empty:
        options = gl_cartesian_product(h1, h2, _alleles_names_)
    else:
        options = gl_cartesian_product_1empty(h1)
    # -------------------- open haplotypes to all the options of allele combinations --------------------/

    # ---------- handle with each allele separately (low to high resolution, ambiguity, etc. ----------/
    full_options = []
    for option in options:
        if not success_gl_string:
            break
        new_option = []
        for i, pair in enumerate(option):
            al1, al2 = pair.split('+')
            success1, al1 = single2high(al1, _alleles_names_[i], idx_fam, idx_par, open_ambiguity_sim,
                                        aux_tools, ard)
            success2, al2 = single2high(al2, _alleles_names_[i], idx_fam, idx_par, open_ambiguity_sim, aux_tools, ard)
            if not success1 or not success2:
                success_gl_string = False
                break
            new_option.append(''.join([al1, '+', al2]))
        new_option = '^'.join(new_option)
        full_options.append(new_option)
    # ---------- handle with each allele separately (low to high resolution, ambiguity, etc. ----------/

    return success_gl_string, full_options


def convert_from_serology(hapF, hapM, alleles_names, group2antigen_dict):
    """
    in serology case, we converted before the antigen to groups (for comparison between children and parents)
    but before we send the gl-strings to GRIMM, we need to reconvert to the antigen (the values that belong to group)
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param alleles_names: alleles_names names
    :param group2antigen_dict: dict that contains the conversion between groups to antigen
    """
    for hap in [hapF.hap1, hapF.hap2, hapM.hap1, hapM.hap2]:
        for allele_name in alleles_names:
            for allele in hap[allele_name]:
                group = ''.join([allele_name, '*', allele])
                # check if the value (group) in dict (low res only. that the reason for condition: if ':' not in allele)
                # if it exists, open to all the serology options and remove the original value
                # (it's only group name, not an antigen)
                # example: A*10 appears, so add A*25, A*25, A*34, A*66 (only the digits) and remove A*10
                if ':' not in allele and group in group2antigen_dict:
                    for ser_value in group2antigen_dict[group]:
                        hap[allele_name].append(ser_value.split('*')[1])


def write_gl_to_file(out_GLstr, gl_string_options, idx_fam, idx_par, races_dict):
    """
    write the gl strings to file
    :param out_GLstr: file path
    :param gl_string_options: the options to gl strings
    :param idx_fam: family index
    :param idx_par: parent index (0: father, 1: mother)
    :param races_dict: dict with data about inserted races to this family
    """
    idx_par = idx_fam + '.0' if idx_par == 0 else idx_fam + '.1'  # if idx_fam=5, so '5.0' to father, '5.1' to mother

    with open(out_GLstr, 'a') as out_file:
        for gl in gl_string_options:
            if 'all' in races_dict:  # Loren simulations. all races are CAU
                out_file.write(''.join([idx_par, ',', gl, ',', 'CAU', ',CAU\n']))
            elif idx_fam in races_dict:  # add the races of the family to the file for GRIMM
                races = ';'.join(races_dict[idx_fam])
                out_file.write(''.join([idx_par, ',', gl, ',', races, ',\n']))
            else:  # no races in this family
                out_file.write(''.join([idx_par, ',', gl, '\n']))


def load_binaries_to_dict(idx_par, idx_fam, single_binary_list, binary_dict):
    idx_par = idx_fam + '.0' if idx_par == 0 else idx_fam + '.1'
    binary_dict[idx_par] = single_binary_list


def create_gl_and_write_to_file(hapF, hapM, alleles_names, idx_fam, par_num, aux_tools, out_GLstr, races_dict,
                                open_ambiguity_sim, errors_in_families):
    """
    after the process in GRAMM, we arrange the options of GL string for each parent and write them to file
    :param hapF: father haplotype
    :param hapM: mother haplotype
    :param alleles_names: alleles names
    :param idx_fam: family index
    :param par_num: number of parents
    :param aux_tools: dict which contains the binaries dict
    :param out_GLstr: file for writing the GL strings
    :param races_dict: dict with data about the races of the families
    :param open_ambiguity_sim: dict with tho options to ambiguity. used in simulations HARD, HARDER, HARDEST
    :param errors_in_families: dict with data about the errors in the families
    """
    binaries_to_GRIMM = []
    ard = pyard.ARD()

    duplicate_hap_if_one_empty(hapF, hapM, aux_tools, idx_fam)

    for idx_par, (h1, h2) in enumerate([(hapF.hap1, hapF.hap2), (hapM.hap1, hapM.hap2)]):
        remove_data_if_just_one_allele_full(h1, h2)
        # create binary list for each parent. the purpose is explained in th function 'create_binary_to_GRIMM'
        binaries_to_GRIMM.append(create_binary_to_GRIMM(h1, h2, alleles_names, par_num))

    for idx_par, (h1, h2) in enumerate([(hapF.hap1, hapF.hap2), (hapM.hap1, hapM.hap2)]):
        success, gl_string_options = create_gl_string(h1, h2, alleles_names, idx_fam, idx_par, par_num,
                                                      open_ambiguity_sim, aux_tools, ard)
        if success:
            write_gl_to_file(out_GLstr, gl_string_options, idx_fam, idx_par, races_dict)
            load_binaries_to_dict(idx_par, idx_fam, binaries_to_GRIMM[idx_par], aux_tools['binary_dict'])
        else:
            errors_in_families[idx_fam] = ['All', '8']
