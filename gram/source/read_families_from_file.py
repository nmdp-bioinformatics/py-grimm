import re
import csv
import json

from .general_aux_funs import convert_to_serology


def get_families(user_file, aux_tools, dont_open_ser=False):
    """
    get families data from user file
    :param user_file: user file
    :param aux_tools: tools for special scenarios. it contains:
        'amb' dict, that keep the ambiguity data, if exists.
        'is_serology' flag, to convert the allele to serology interpretation, if the user choose that
        'ser_dict' that contains the serology conversion rules
    :param dont_open_ser: True when we read the families in "create_results", so we do not want to convert to serology
    (because we do it into the code, in different way (only children)
    :return: dict with data about the families, in this structure:
        keys: families indexes (1, 2..)
        values: dict for each family
            in these sub-dicts (for each family):
                keys: indexes of family members: F, M, 1, 2 ...
                values: alleles_names data
                    in these sub-sub dict (of alleles_names data for each family member):
                    keys: alleles_names names (A, B..)
                    values: list of numbers of each alleles_names (e.g: [02:01, 30:04])
    """
    file = open(user_file, 'r')
    reader = csv.reader(file)
    next(reader, None)  # skip headers
    families_dict = {}

    for line in reader:
        if not line:  # empty line
            continue
        id_family = line[0]
        id_person = line[1]
        indv_alleles, data_exist = get_individual(line, aux_tools, id_family, dont_open_ser)
        if data_exist:
            if id_family not in families_dict:
                families_dict[id_family] = {}  # create a sub dict to the family in families_dict
            families_dict[id_family][id_person] = indv_alleles

    file.close()

    return families_dict


def get_individual(line, aux_tools, id_family, dont_open_ser):
    """
    get individual data
    :param line: line in file
    :param aux_tools: explained in "get_families" function
    :param id_family: family index
    :param dont_open_ser: True when we read the families in "create_results", so we do not want to convert to serology
    :return: individual data, and flag that equal to False if there is no data about this individual
    """
    indv_alleles = {}
    data_exist = True

    alleles_map = {2: 'A', 4: 'B', 6: 'C', 8: 'DRB1', 10: 'DQB1'}
    for i in range(2, 11, 2):
        alleles_pair = line[i:i+2]
        new_alleles_pair = process_alleles(alleles_pair, alleles_map[i], aux_tools, id_family, dont_open_ser)
        indv_alleles[alleles_map[i]] = new_alleles_pair
    if not any(item for sublist in indv_alleles.values() for item in sublist):  # completely empty: {('',''),('','')..}
        data_exist = False

    return indv_alleles, data_exist


def process_alleles(pair, allele, aux_tools, id_family, dont_open_ser):
    """
    process alleles_names: replace irrelevant characters, remove ambiguity (if exists) and save it to a dict,
    convert to serology if needed
    :param pair: alleles_names pair
    :param allele: the specific allele of this pair (A/B..)
    :param aux_tools: explained in "get_families" function
    :param id_family: family index
    :param dont_open_ser: True when we read the families in "create_results", so we do not want to convert to serology
    :return: alleles_names pair, after process
    """
    amb = aux_tools['amb']  # amb is a dict of ambiguity of alleles_names
    is_serology = aux_tools['is_serology']  # is_serology is boolean (serology samples or genetic)
    ser_dict = aux_tools['antigen2group']  # ser_dict is a dict for the conversion to serology interpretation

    al1, al2 = pair[0], pair[1]

    # in simulation files, al1='02:01', al2='', mean that it homozygous (they both 02:01)
    if al2 == '':
        al2 = al1

    new_als = []
    for al in [al1, al2]:
        al = str(al).replace('p', '')  # is necessary? after we do re.sub to [a-zA-Z]

        # more than 4 digits: 01:02:01 -> 01:02. we don't remove from the sixth char, because could be letters: 01:APCJ
        if al.count(':') > 1:
            parts = al.split(':')
            al = ''.join([parts[0], ':', parts[1]])

        """
            Save the ambiguity in dict and remove it from al1, al2.
            That because in the stage of comparing between parents and children, it could disturb.
            For example: parent: A*02:BJFV, _child_:A*02:02, we will not recognize it could match.
            Before the insertion to GRIMM, we restore it to data.
        """
        if ":" in al and al.split(':')[1].isupper():
            if id_family not in amb:
                amb[id_family] = {}
            amb[id_family][allele + '*' + al.split(':')[0]] = al.split(':')[1]

        al = re.sub('[a-zA-Z]', '', al)  # remove ambiguity

        if len(al) == 1 and al.isdigit():  # one digit to two digits (2 -> 02)
            al = '0' + al

        if is_serology and not dont_open_ser:
            al = convert_to_serology(ser_dict, allele, al)

        new_als.append(al)

    return new_als

