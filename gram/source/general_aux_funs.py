

def equal_als(al_1, al_2):
    """
    :param al_1: first allele
    :param al_2: second allele
    :return: true if the allele are identical or including, false otherwise
    """
    if al_1 == al_2 or al_1 is al_2 or al_1.startswith(al_2) or al_2.startswith(al_1):
        if bool(al_1) == bool(al_2):
            return True
    return False


def empty_hap(hap):
    """
    check if haplotype is empty (no data, just: {'A':[], 'B':[], ..}
    @param hap: dict represents haplotype
    @return: True if empty, else False
    """
    for value in hap.values():
        if value:
            return False
    return True


def convert_to_serology(ser_dict, allele_name, single_al):
    """
    when data is serology, the alleles_names need to be converted, for example:
    ser_dict = {"A*23": "A*09", "A*24": "A*09"..} so if we get allele = A, single_al = 23, change to 09
    :param ser_dict: serology dict
    :param allele_name: allele name
    :param single_al: single allele
    :return:
    """
    # case of low res (A*23 -> A*09)
    if allele_name + '*' + single_al in ser_dict:
        return ser_dict[allele_name + '*' + single_al].split('*')[1]
    # case of high res (A*23:01 -> A*09)
    if ':' in single_al and allele_name + '*' + single_al.split(':')[0] in ser_dict:
        return ser_dict[allele_name + '*' + single_al.split(':')[0]].split('*')[1]
    return single_al
