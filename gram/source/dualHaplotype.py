from .als_update import Als


class DualHaplotype:
    """
    this class represent 2 haplotypes of a parent (father or mother)
    each haplotype is a dictionary, which contains alleles_names data
    """
    def __init__(self, alleles_names):
        """
        create dict for each hoplotype, with empty Als (=special list for alleles_names) for each allele: {A: [], B: [], ...}
        :param alleles_names: alleles_names
        """
        self.hap1 = {}
        self.hap2 = {}
        for al_name in alleles_names:
            self.hap1[al_name] = Als()
            self.hap2[al_name] = Als()

    def insert_parents_data(self, family, parent, par_num):
        """
        insertion alleles_names data to parents haplotypes
        first insertion - permanent (hap1: {A:[02]...}, hap2: {A:[03]...})
        others - two options (hap1: {A:[02], B:[01, 08]...}, hap2: {A:[03], B:[01, 08]...})
        (if there is an homozygous allele - it's permanent, too)
        :param family: dict with family data
        :param parent: F or M
        :param par_num: parents number: 0/1/2
        """
        is_permanent = True  # after one permanent insertion, change to False

        for allele_name in family[parent]:
            al1 = family[parent][allele_name][0]
            al2 = family[parent][allele_name][1]
            if al1 == al2:  # homozygous
                self.hap1[allele_name].append(al1)
                self.hap2[allele_name].append(al2)
            else:
                # if there wasn't permanent insertion yet,
                # and we are not in a situation that there are 2 parents with same alleles_names (it will be problem in the
                # stage that we associate children with parents haplotypes)
                # -> so we do permanent insertion
                if is_permanent and not(par_num == 2 and family['F'][allele_name] == family['M'][allele_name]):
                    self.hap1[allele_name].append(al1)
                    self.hap2[allele_name].append(al2)
                    is_permanent = False
                else:
                    self.hap1[allele_name].extend([al1, al2])
                    self.hap2[allele_name].extend([al1, al2])
