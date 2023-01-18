import sys
import os
import math

# import matplotlib.pyplot as plt
from .em_util import *

sys.path.insert(0, os.path.join(".."))
from grim.imputation.impute import Imputation
from grim.imputation.impute import clean_up_gl


class algo:
    def __init__(self, config, pop):
        self.freq_file = config["freq_file"]
        self.f_data = config["imputation_input_file"]
        self.imputation_res = config["imputation_out_hap_freq_file"]
        self.pop = pop
        self.loci_map = config["loci_map"]
        self.all_loci = list(config["loci_map"].keys())
        self.haps_num = []
        self.log_muugs = []
        self.freq_change = []
        self.max_size = config["memory_max"]
        self.min_size = config["memory_min"]
        self.cutoff_init = config["cutoff_init"]

    ## count number appearance of all haplotypes, and number appearance of each haplotype
    def add_haps_to_dict(self, dict_loci, phases_list):
        sum = 0
        for phases in phases_list:
            for j in phases:
                for i in j:
                    for k in i:
                        phase = ("~").join(k)
                        if not phase in dict_loci:
                            # if len(dict_loci) > self.max_size:
                            # dict_loci = cut_dict(dict_loci, self.min_size)
                            dict_loci[phase] = 1
                        else:
                            dict_loci[phase] += 1
                        sum += 1
        return sum  # , dict_loci

    ## calculate the probability of each haplotype and write to file - hap,pop,prob
    def write_freqs_to_file(self, dict_loci, sum):
        f_freq = open(self.freq_file, "w")
        for hap in dict_loci:
            f_freq.write(hap + "," + self.pop + "," + str(dict_loci[hap] / sum) + "\n")

        f_freq.close()

    def create_guess(self):

        haps_count = 0
        imputation = Imputation()
        dict_loci = {}
        num_samples = 0

        with open(self.f_data) as lines:
            for name_gl in lines.readlines():
                if "," in name_gl:
                    name_gl = name_gl.split(",")
                else:
                    name_gl = name_gl.split("%")
                name_gl[1] = clean_up_gl(name_gl[1])
                # name_gl = changeFormat(name_gl)
                num_samples += 1
                if self.is_full_haplo(name_gl[1]):
                    if len(dict_loci) > self.max_size:
                        dict_loci = cut_dict(dict_loci, self.min_size)

                    # print(str(num_samples) + ' ' + str(len(dict_loci)))
                    phases_list = imputation.open_gl_string(
                        name_gl[1], self.cutoff_init
                    )
                    if not phases_list:
                        continue
                    # sum, dict_loci = self.add_haps_to_dict(dict_loci, phases_list)
                    # haps_count += sum
                    haps_count += self.add_haps_to_dict(dict_loci, phases_list)

        lines.close()

        self.haps_num.append(len(dict_loci))

        # calculate the probability of each haplotype and write to file
        self.write_freqs_to_file(dict_loci, haps_count)

        # update epsilon
        self.eps_conv = 1 / (num_samples)
        # self.eps_conv = 0.00001
        return self.eps_conv, len(dict_loci), num_samples

    ## normalize the genotypes sum of each person to 1
    ## calc Hi - sum(G(i,j) +2G(i,i)) * 0.5
    ## normalize the all haplotypes sum to 1
    def calc_new_prob(self):
        dict_probs = {}
        sum_all_geno = 0
        # imputation_res = 'output/' + self.f_data + '_outMPG.frq'

        j = 0
        sum_log = 0
        with open(self.imputation_res) as imputat_res:
            res = imputat_res.readline()
            while res:
                id = res.split(",")[0]
                sum = 0
                list_person = []

                # find the all results of single person
                # sum the prob of single person
                while res and res.split(",")[0] == id:
                    res = res.strip().split(",")
                    p = float(res[2])
                    hap1, hap2 = res[1].split("+")
                    list_person.append([hap1, hap2, p])
                    sum += p
                    res = imputat_res.readline()

                sum_log += math.log(sum)
                if len(dict_probs) > self.max_size:
                    dict_probs = cut_dict(dict_probs, self.min_size)
                # add to each haplotype's prob the normalize prob from each person
                # sum the total prob of all persons
                for haps in list_person:
                    p = haps[2] / sum  # normalize prob by person
                    sum_all_geno += p  # sum the total prob
                    p *= 0.5
                    if haps[0] in dict_probs:
                        dict_probs[haps[0]] += p
                    else:
                        dict_probs[haps[0]] = p

                    if haps[1] in dict_probs:
                        dict_probs[haps[1]] += p
                    else:
                        dict_probs[haps[1]] = p

        print(sum_all_geno)

        self.haps_num.append(len(dict_probs))
        self.log_muugs.append(sum_log)

        self.write_freqs_to_file(dict_probs, sum_all_geno)

        return sum_log, len(dict_probs)

    # if all the difference
    def check_converges(self, dict_old):
        not_conv = False
        sum_freq_change = 0.0
        print(self.eps_conv)
        with open(self.freq_file) as new_freq_file:
            hap = new_freq_file.readline()
            while hap:
                hap, pop, prob = hap.strip().split(",")
                prob = float(prob)
                if hap in dict_old:
                    sum_freq_change += (prob - dict_old[hap]) ** 2
                    if abs(prob - dict_old[hap]) > self.eps_conv:
                        # return True
                        not_conv = True
                elif prob > self.eps_conv:
                    # return True
                    not_conv = True
                hap = new_freq_file.readline()

        self.freq_change.append(sum_freq_change)
        return not_conv, sum_freq_change
        # return False

    """def plots(self):
        Plot(self.log_muugs[1:], 'Log likelihood - EM1', 'sum(log(sum(MUUGS)))', 'Iteration', self.plot_path + 'logEM1.png' )
        Plot(self.haps_num[1:], 'Number of haplotypes per iteration - EM1', 'Haplotypes number', 'Iteration', self.plot_path + 'hapnumEM1.png' )
        Plot(self.freq_change, 'Change in frequencies - EM1', ' sum((P(Hi) - P(Hi-1))^2)', 'Iteration', self.plot_path + 'freqEM1.png' )"""

    def calc_log_likelihood(self):
        dict_probs = {}
        sum_all_geno = 0
        # imputation_res = 'output/' + self.f_data + '_outMPG.frq'

        j = 0
        sum_log = 0
        with open(self.imputation_res) as imputat_res:
            res = imputat_res.readline()
            while res:
                id = res.split(",")[0]
                sum = 0
                list_person = []

                # find the all results of single person
                # sum the prob of single person
                while res and res.split(",")[0] == id:
                    res = res.split(",")
                    p = float(res[3])
                    list_person.append([res[1], res[2], p])
                    sum += p
                    res = imputat_res.readline()

                sum_log += math.log(sum)
        return sum_log

    """def is_full_haplo(self, hap):
        if all(loci in hap for loci in self.all_loci):
            return True
        return False"""

    def is_full_haplo(self, hap):
        list_loci_in_hap = []
        all_loci = hap.split("^")
        for locus in all_loci:
            list_loci_in_hap.append(self.loci_map[locus.split("*")[0]])
        if all(loci in list_loci_in_hap for loci in list(self.loci_map.values())):
            return True
        return False
