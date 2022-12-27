CONFIG = {
    "pops": [
        "AAFA", "AFB", "AINDI", "AISC", "ALANAM", "AMIND", "CARB", "CARHIS", "CARIBI",
        "HAWI", "FILII", "KORI", "JAPI", "MSWHIS", "MENAFC", "NAMER", "NCHI", "SCAHIS",
        "SCAMB", "SCSEAI", "VIET", "AFA", "API", "CAU", "HIS", "NAM", "Kavkaz", "Ashkenazi",
        "Ethiopian", "Sephardi", "Others", "Arab", "all"
    ],
    "freq_trim_threshold": 0.5,
    "priority": {
        "alpha": 0.4999999,
        "eta": 0,
        "beta": 1e-7,
        "gamma": 1e-7,
        "delta": 0.4999999
    },
    "UNK_priors": "SR",
    "FULL_LOCI": "ABCQR",
    "loci_map": {
        "A": 1,
        "B": 2,
        "C": 3,
        "DQB1": 4,
        "DRB1": 5
    },
    "factor_missing_data": 0.0001,
    "matrix_planb": [
        [[1, 2, 3, 4, 5]],
        [[1, 2, 3], [4, 5]],
        [[1], [2, 3], [4, 5]],
        [[1, 2, 3], [4], [5]],
        [[1], [2, 3], [4], [5]],
        [[1], [2], [3], [4], [5]]
    ],
    "planb": True,
    "number_of_options_threshold": 100000,
    "epsilon": 1e-3,
    "number_of_results": 10,
    "number_of_pop_results": 100,
    "output_MUUG": True,
    "output_haplotypes": False,
    "max_haplotypes_number_in_phase": 300,
    "bin_imputation_input_file": "None",
    "use_pops_count_file": False,
    "nodes_for_plan_A": [],
    "save_mode": False

    # fields that must be given by the user:

    # "imputation_in_file": input file for the imputation.
    # "imputation_out_umug_freq_filename": genotype frequencies output file. the only file relevant for the matching.
    # "imputation_out_umug_pops_filename": genotype pops output file
    # "imputation_out_hap_freq_filename": haplotype frequencies output file.
    # "imputation_out_hap_pops_filename": haplotype pops output file.
    # "imputation_out_miss_filename": miss output file.
    # "imputation_out_problem_filename": problem output file.
}


