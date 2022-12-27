from __future__ import annotations

from datetime import datetime
from typing import Iterable

from GRMA.Utilities.cutils import cdrop_less_than_7_matches, ccheck_similarity

from collections.abc import Sequence

LENGTH_OF_NUMBERS_IN_ALLELE: int = 4
TO_SEROLOGY: dict[str, int] = {"A2": 1, "A203": 1, "A210": 1, "A9": 2, "A23": 2, "A24": 5, "A10": 3, "A25": 3, "A26": 3,
                               "A34": 3,
                               "A66": 3, "A19": 4,
                               "A29": 4, "A30": 4, "A31": 4, "A32": 4, "A33": 4, "A74": 4, "A2403": 5, "A28": 6,
                               "A68": 6, "A69": 6,
                               "B5": 7,
                               "B51": 19, "B52": 7, "B7": 8, "B703": 8, "B12": 9, "B44": 9, "B45": 9, "B14": 10,
                               "B64": 10,
                               "B65": 10, "B15": 11,
                               "B62": 11, "B63": 11, "B75": 11, "B76": 11, "B77": 11, "B16": 12, "B38": 12, "B39": 17,
                               "B17": 13,
                               "B57": 13,
                               "B58": 13, "B21": 14, "B49": 14, "B50": 14, "B4005": 14, "B22": 15, "B54": 15, "B55": 15,
                               "B56": 15,
                               "B27": 16,
                               "B2708": 16, "B3901": 17, "B3902": 17, "B40": 18, "B60": 18, "B61": 18, "B5102": 19,
                               "B5103": 19,
                               "B70": 20,
                               "B71": 20, "B72": 20, "C3": 21, "C9": 21, "C10": 21, "DQ1": 22, "DQ5": 22, "DQ6": 22,
                               "DQ3": 23,
                               "DQ7": 23, "DQ8": 23,
                               "DQ9": 23, "DR1": 24, "DR103": 24, "DR2": 25, "DR15": 25, "DR16": 25, "DR3": 26,
                               "DR17": 26,
                               "DR18": 26, "DR5": 27,
                               "DR11": 27, "DR12": 27, "DR6": 28, "DR13": 28, "DR14": 29, "DR1403": 29, "DR1404": 29}


def list_to_genotype(geno_list: Sequence[int]) -> str:
    """convert genotype from list[int] to str according to gl-string format"""
    allele_dict = {
        0: "A*", 1: "A*",
        2: "B*", 3: "B*",
        4: "C*", 5: "C*",
        6: "DQB1*", 7: "DQB1*",
        8: "DRB1*", 9: "DRB1*"
    }
    geno: str = ""
    for i in range(10):
        allele_num = str(geno_list[i])
        geno += f"{allele_dict[i]}"
        geno += f"{allele_num[:-2]}:{allele_num[-2:]}".zfill(5)
        geno += "+" if i % 2 == 0 else "^"
    return geno[:-1]


def donor_mismatch_format(don_geno: Sequence[int], pat_geno: Sequence[int]) -> str:
    """
    Takes patient's genotype and donor's genotype. \n
    Replace donor's alleles places to match to patient's genotypes. \n
    If there is a mismatch:
     * surround allele with [] for match in serology.
     * surround allele with () for allele mismatch.
    Returns the donor's genotype with its alleles separated by commas.
    """
    global TO_SEROLOGY, LENGTH_OF_NUMBERS_IN_ALLELE
    allele_dict = {
        0: "A", 1: "A",
        2: "B", 3: "B",
        4: "C", 5: "C",
        6: "DQ", 7: "DQ",
        8: "DR", 9: "DR"
    }
    reformatted_donor_geno = ['' for _ in range(10)]

    def allele_style(patients_alleles: int, match: bool,
                     alleles_place: int, donors_alleles: int | None = None) -> str:
        """
        Replace donor's alleles places to match to patient's genotypes. \n
        If there is a mismatch:
         * surround allele with [] for match in serology.
         * surround allele with () for allele mismatch.
        """
        patients_allele_to_str = [str(patients_alleles).zfill(LENGTH_OF_NUMBERS_IN_ALLELE)[:2],
                                  str(patients_alleles).zfill(LENGTH_OF_NUMBERS_IN_ALLELE)[2:]]

        donors_allele_to_str = [str(donors_alleles).zfill(LENGTH_OF_NUMBERS_IN_ALLELE)[:2],
                                str(donors_alleles).zfill(LENGTH_OF_NUMBERS_IN_ALLELE)[
                                2:]] if donors_alleles is not None else None
        if match:
            return f"{patients_allele_to_str[0]}:{patients_allele_to_str[1]}"

        if TO_SEROLOGY.get(f"{allele_dict[alleles_place]}{patients_allele_to_str[0]}",
                           f"{allele_dict[alleles_place]}{patients_allele_to_str[0]}") == \
                TO_SEROLOGY.get(f"{allele_dict[alleles_place]}{donors_allele_to_str[0]}",
                                f"{allele_dict[alleles_place]}{donors_allele_to_str[0]}"):
            return f"[{patients_allele_to_str[0]}:{patients_allele_to_str[1]}]"

        else:
            return f"({patients_allele_to_str[0]}:{patients_allele_to_str[1]})"

    for i in range(0, 10, 2):
        if don_geno[i] == pat_geno[i]:
            reformatted_donor_geno[i] = allele_style(don_geno[i], True, i)
            reformatted_donor_geno[i + 1] = allele_style(don_geno[i + 1], don_geno[i + 1] == pat_geno[i + 1], i,
                                                         pat_geno[i + 1])

        elif don_geno[i + 1] == pat_geno[i + 1]:
            reformatted_donor_geno[i] = allele_style(don_geno[i], False, i, pat_geno[i])
            reformatted_donor_geno[i + 1] = allele_style(don_geno[i + 1], True, i)

        elif don_geno[i] == pat_geno[i + 1]:
            reformatted_donor_geno[i] = allele_style(don_geno[i + 1], don_geno[i + 1] == pat_geno[i], i,
                                                     pat_geno[i])
            reformatted_donor_geno[i + 1] = allele_style(don_geno[i], True, i)

        elif don_geno[i + 1] == pat_geno[i]:
            reformatted_donor_geno[i] = allele_style(don_geno[i + 1], True, i)
            reformatted_donor_geno[i + 1] = allele_style(don_geno[i], False, i, pat_geno[i + 1])

        else:
            reformatted_donor_geno[i] = allele_style(don_geno[i], False, i, pat_geno[i])
            reformatted_donor_geno[i + 1] = allele_style(don_geno[i + 1], False, i, pat_geno[i + 1])

    return ",".join(reformatted_donor_geno)


def tuple_geno_to_int(tuple_geno: Iterable[int]):
    """convert genotype from Iterable to integer"""
    string = ""
    fill: int = 4
    for i in tuple_geno:
        string += str(i).zfill(fill)
    return int(string)


def print_time(log: str):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(f"[{current_time}] -- {log}")


def drop_less_than_7_matches(ids, similarities):
    return cdrop_less_than_7_matches(ids, similarities)


def check_similarity(patients_geno, donors_genos, allele_range, init_count_similar):
    return ccheck_similarity(patients_geno, donors_genos, allele_range, init_count_similar)


def gl_string_to_integers(genotype: str) -> Sequence[int]:
    genotype = genotype.replace('+', '~').replace('^', '~').replace(":", "")
    genotype = [int(allele.split("*")[1][:4]) for allele in genotype.split("~")]
    return genotype
