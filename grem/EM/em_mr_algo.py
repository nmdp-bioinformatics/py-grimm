import numpy as np
from .em_util import *
import math


class algo_mr:
    def __init__(self, config, epsilon):

        self.freq_file = config["freq_file"]
        self.f_data = config["imputation_input_file"]
        self.imputation_res = config["imputation_out_hap_freq_file"]
        self.pops = config["pops"]
        self.eps_conv = epsilon

        self.haps_num = []
        self.log_muugs = []
        self.freq_change = []

        self.max_size = config["memory_max"]
        self.min_size = config["memory_min"]

    ## count number appearance of all haplotypes, and number appearance of each haplotype
    def count_haps_by_race(self, dict_haps, phases_list, race1, race2):

        for phases in phases_list:
            for j in phases:
                for i in j:
                    for k in i:
                        phase = ("~").join(k)
                        if phase in dict_haps:
                            dict_haps[phase][0][self.pops.index(race1)] += 0.5
                            dict_haps[phase][0][self.pops.index(race2)] += 0.5

    # normalize prob for each race and write to file
    def write_freqs_guess_to_file(self, dict_haps):
        count_list_by_prob = np.zeros(len(self.pops))

        sum_all = 0
        f_freq = open(self.freq_file, "w")
        for hap, probs in dict_haps.items():
            sum = np.sum(probs[0])
            if sum != 0:
                prob = probs[0] / sum * float(probs[1])
                for pop in self.pops:
                    p = prob[self.pops.index(pop)]
                    if p != 0:
                        f_freq.write(hap + "," + pop + "," + str(p) + "\n")
                        count_list_by_prob[self.pops.index(pop)] += p

        print("sum_all")
        print(np.sum(count_list_by_prob))
        f_freq.close()
        return count_list_by_prob

    def create_guess(self):

        freqs_list = []
        with open(self.freq_file) as hpf_file:
            for line in hpf_file:

                for pop in self.pops:
                    l = line.strip().split(",")
                    l[1] = pop
                    # l[2] = str(float(l[2]))
                    freqs_list.append((",").join(l) + "\n")

        file_freq_out = open(self.freq_file, "w")
        for freq in freqs_list:
            file_freq_out.write(freq)

    ## calculate the probability of each haplotype and write to file - hap,pop,prob
    # normalize each pop prob to 1

    def write_freqs_to_file(self, dict_loci, sum):
        f_freq = open(self.freq_file, "w")
        for hap, probs in dict_loci.items():
            for pop in self.pops:
                p = probs[0][self.pops.index(pop)] / sum[self.pops.index(pop)]
                if p != 0:
                    f_freq.write(hap + "," + pop + "," + str(p) + "\n")

        f_freq.close()

    ## normalize the genotypes sum of each person to 1
    ## calc Hi - sum(G(i,j) +2G(i,i)) * 0.5
    ## normalize the all haplotypes sum to 1
    def calc_new_prob(self):
        dict_probs = {}
        sum_all_geno = 0
        pops_len = len(self.pops)

        j = 0
        sum_log = 0
        count_list_by_pop = np.zeros(len(self.pops))
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

                if len(dict_probs) > self.max_size:
                    dict_probs = self.cut_dict(dict_probs)

                # add to each haplotype's prob the normalize prob from each person
                # sum the total prob of all persons
                for haps in list_person:
                    p = haps[2] / sum  # normalize prob by person
                    sum_all_geno += p  # sum the total prob
                    p *= 0.5
                    for i in range(2):
                        hap, race = haps[i].split(";")
                        if hap in dict_probs:
                            dict_probs[hap][1] += p
                        else:
                            dict_probs[hap] = [np.zeros(pops_len), p]
                        dict_probs[hap][0][self.pops.index(race)] += p

                    """if haps[1] in dict_probs:
                        dict_probs[haps[1]] += p
                    else:
                        dict_probs[haps[1]] = p"""

        # normalize the all haplotypes sum to 1 and sum prob in each population
        for hap, probs in dict_probs.items():
            for pop in self.pops:
                probs[0][self.pops.index(pop)] = (
                    probs[0][self.pops.index(pop)] / sum_all_geno
                )
                count_list_by_pop[self.pops.index(pop)] += probs[0][
                    self.pops.index(pop)
                ]  ##/sum_all_geno

        print(sum_all_geno)
        sum_log = sum_log / self.num_sample
        self.haps_num.append(len(dict_probs))
        self.log_muugs.append(sum_log)

        self.write_freqs_to_file(dict_probs, count_list_by_pop)

        return count_list_by_pop, sum_log

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
                    sum_freq_change += (prob - dict_old[hap][self.pops.index(pop)]) ** 2
                    if abs(prob - dict_old[hap][self.pops.index(pop)]) > self.eps_conv:
                        # return True
                        not_conv = True
                elif prob > self.eps_conv:
                    sum_freq_change += prob**2
                    # return True
                    not_conv = True
                hap = new_freq_file.readline()

        self.freq_change.append(sum_freq_change)
        return not_conv

    """def plots(self):
        Plot(self.log_muugs[1:], 'Log likelihood - EM2', 'sum(log(sum(MUUGS)))', 'Iteration', self.plot_path + 'logEM2.png' )
        Plot(self.haps_num[1:], 'Number of haplotypes per iteration - EM2', 'Haplotypes number', 'Iteration', self.plot_path + 'hapnumEM2.png' )
        Plot(self.freq_change, 'Change in frequencies - EM2', ' sum((P(Hi) - P(Hi-1))^2)', 'Iteration', self.plot_path + 'freqEM2.png' )"""

    def cut_dict(self, dict_val):
        cut_dict = {}
        dict_val = sorted(dict_val.items(), key=lambda x: x[1][1], reverse=True)
        # dict_val.clear()
        while len(cut_dict) < self.min_size:
            cut_dict[dict_val[0][0]] = dict_val[0][1]
            del dict_val[0]

        return cut_dict

    def create_eps(self):

        haps_count = 0
        # imputation = Imputation()
        dict_loci = {}
        num_samples = 0

        with open(self.f_data) as lines:
            for name_gl in lines.readlines():

                num_samples += 1
        lines.close()

        # update epsilon
        self.eps_conv = 1 / (num_samples)
        self.num_sample = num_samples
        return self.eps_conv

    def calc_log_likelihood(self):
        dict_probs = {}
        sum_all_geno = 0
        pops_len = len(self.pops)

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


# D009116 - 1.0480316388972355e-06
#'D008326' - 7.8677600616405e-14
#'D000605' - 4.616528288275305e-06
