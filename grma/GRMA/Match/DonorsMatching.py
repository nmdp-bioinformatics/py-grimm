from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import List, Tuple, Set, Iterable, Dict

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from GRMA.Match.GraphWrapper import Graph
from GRMA.Utilities.geno_representation import HashableArray, ClassMinusOne
from GRMA.Utilities.utils import donor_mismatch_format, \
    drop_less_than_7_matches, check_similarity, gl_string_to_integers, tuple_geno_to_int

DONORS_DB: pd.DataFrame = pd.DataFrame()
ZEROS: HashableArray = HashableArray([0])
ALLELES_IN_CLASS_I: int = 6
ALLELES_IN_CLASS_II: int = 4


def set_Database(donors_db: pd.DataFrame = pd.DataFrame()):
    """
    Set a database for a search.
    Use this function before the matching if you wish to add fields for the result df.
    """
    global DONORS_DB
    DONORS_DB = donors_db


def _init_results_df(donors_info):
    """Initialize matching donors' df"""
    global DONORS_DB
    fields_in_results = {
        "Patient_ID": [], "Donor_ID": [],
        "Number_Of_Mismatches": [], "Matching_Probability": [],
        "Match_Probability_A_1": [], "Match_Probability_A_2": [],
        "Match_Probability_B_1": [], "Match_Probability_B_2": [],
        "Match_Probability_C_1": [], "Match_Probability_C_2": [],
        "Match_Probability_DQB1_1": [], "Match_Probability_DQB1_2": [],
        "Match_Probability_DRB1_1": [], "Match_Probability_DRB1_2": [],
        "Permissive/Non-Permissive": []
    }

    donors_db_fields = DONORS_DB.columns.values.tolist()
    for di in donors_info:
        if di in donors_db_fields:
            fields_in_results[di] = []
    return pd.DataFrame(fields_in_results)


class DonorsMatching(object):
    """DonorsMatching class is in charge of the matching process"""
    __slots__ = "_graph", "_patients_graph", "patients"

    def __init__(self, graph: Graph):
        self._graph: Graph = graph
        self._patients_graph: nx.DiGraph = nx.DiGraph()
        self.patients: dict[int, Sequence[int]] = {}

    # <print matching information to csv>
    def print_most_common_genotype(self, don_id: int, pat_geno: Sequence[int]) -> str:
        """Takes a donor ID and a genotype. \n
        Returns the mismatch format of the most common genotype of the donor."""
        don_geno = ""
        geno_max_prob = 0
        for geno in self._graph.neighbors(don_id):
            if geno[1] > geno_max_prob:
                geno_max_prob = geno[1]
                don_geno = geno[0]
        return donor_mismatch_format(don_geno, pat_geno)

    def probability_to_allele(self, don_id: int, pat_geno: Sequence[int]) -> List[float]:
        """Takes a donor ID and a genotype. \n
        Returns the probability of match for each allele"""
        probs = [0 for _ in range(10)]

        for i, allele in enumerate(pat_geno):
            p = 0
            for don_geno, don_weight in self._graph.neighbors(don_id):
                if allele in don_geno:
                    p += don_weight
            probs[i] = int(round(p * 100))

        return probs

    # </print matching information to csv>

    def find_genotype_from_subclass(self, sub: int) -> np.ndarray:
        """Takes an integer subclass. \n
        Returns the genotypes which are connected to it in the graph"""
        return self._graph.neighbors_2nd(sub)

    def find_donor_from_geno(self, geno_id: int) -> Sequence[int]:
        """Get the LOL ID of a genotype. \n
        Return its neighbors - all the donors that has this genotype."""
        donors = []
        for iid, w in self._graph.neighbors(geno_id, search_lol_id=True):
            donors.append(iid)
        return donors

    def find_candidates(self, subclass: ClassMinusOne, genos: Iterator) -> None:
        """
        Takes a subclass object and an iterator of patients' genotypes that has this subclass. \n
        Add to patients_graph all the genotypes candidates.
        """
        # check only the locuses that are not certain to match
        if subclass.class_num == 0:
            allele_range_to_check = np.array([6, 8, subclass.allele_num], dtype=np.uint8)
        else:
            allele_range_to_check = np.array([0, 2, 4, subclass.allele_num], dtype=np.uint8)

        # number of alleles that already match due to match in subclass
        matched_alleles: int = (ALLELES_IN_CLASS_II if subclass.class_num == 1 else ALLELES_IN_CLASS_I) - 2
        # get all candidates
        genotypes_id_from_subclass, genotypes_value_from_subclass = self.find_genotype_from_subclass(subclass.subclass)

        for geno in genos:
            # check similarity between geno and all the candidates
            similarities = check_similarity(geno.np(),
                                            genotypes_value_from_subclass, allele_range_to_check,
                                            matched_alleles)
            # remove candidates with less than 7 matches
            candidates_to_iterate = drop_less_than_7_matches(genotypes_id_from_subclass, similarities)
            for candidate_id, similarity in candidates_to_iterate:
                # iterate over all the patients with the genotype
                for patient_id in self._patients_graph.neighbors(geno):
                    # patient's geno index
                    geno_num = self._patients_graph[geno][patient_id]["geno_num"]
                    # patient's geno probability
                    probability = self._patients_graph[geno][patient_id]["probability"]

                    # add the genotype id as a neighbor to the patient
                    if candidate_id in self._patients_graph.adj[patient_id]:
                        self._patients_graph[patient_id][candidate_id]['weight'][geno_num] = [probability,
                                                                                              similarity]
                    else:
                        if candidate_id not in self._patients_graph:
                            self._patients_graph.add_node(candidate_id)
                        self._patients_graph.add_edge(patient_id, candidate_id,
                                                      weight={geno_num: [probability, similarity]})

    def add_subclasses(self, genotype: HashableArray, subclasses: list[ClassMinusOne]) -> list[ClassMinusOne]:
        """Takes a genotype, and a list of subclasses that were already processed,
        process the subclasses of the given genotype, append them to the given list and returns the list."""
        classes = [genotype[:ALLELES_IN_CLASS_I], genotype[ALLELES_IN_CLASS_I:]]
        num_of_alleles_in_class = [ALLELES_IN_CLASS_I, ALLELES_IN_CLASS_II]

        # class one is considered as 0.
        # class two is considered as 1.
        CLASSES = [0, 1]
        for class_num in CLASSES:
            for k in range(0, num_of_alleles_in_class[class_num]):
                # set the missing allele to always be the second allele in the locus
                if k % 2 == 0:
                    sub = tuple_geno_to_int(classes[class_num][0: k] + ZEROS + classes[class_num][k + 1:])
                else:
                    sub = tuple_geno_to_int(classes[class_num][0: k - 1] + ZEROS +
                                            classes[class_num][k - 1: k] + classes[class_num][k + 1:])

                # missing allele number is the index of the first allele of the locus the missing allele belongs to.
                # Could be [0, 2, 4, 6, 8]
                missing_allele_num = ALLELES_IN_CLASS_I * class_num + 2 * (k // 2)
                subclass = ClassMinusOne(subclass=sub,
                                         class_num=class_num,
                                         allele_num=missing_allele_num)

                # add subclass -> genotype edge to patients graph
                if subclass not in self._patients_graph:
                    subclasses.append(subclass)
                self._patients_graph.add_edge(subclass, genotype)
        return subclasses

    def create_patients_graph(self, f_patients: str) -> list[ClassMinusOne]:
        """
        create patients graph. \n
        *takes in consideration that grimm outputs for each patient different genotypes*
        """
        self._patients_graph: nx.DiGraph = nx.DiGraph()
        prob_dict: dict = {}  # {geno: [i, prob]}
        total_prob: float = 0
        last_patient: int = -1
        subclasses: list[ClassMinusOne] = []

        for line in open(f_patients).readlines():
            # retrieve all line's parameters
            patient_id, geno, prob, index = line.strip().split(',')

            geno = gl_string_to_integers(geno)
            patient_id = int(patient_id)
            index = int(index)
            prob = float(prob)

            # handle new patient appearance in file
            if index == 0:
                # set normalized probabilities
                for HLA, (_, probability) in prob_dict.items():
                    self._patients_graph.edges[HLA, last_patient]['probability'] = \
                        probability / total_prob

                # initialize parameters
                prob_dict = {}
                total_prob = 0
                self.patients[patient_id] = geno
                last_patient = patient_id

            # sort alleles for each HLA-X
            for x in range(0, 10, 2):
                geno[x: x + 2] = sorted(geno[x: x + 2])
            geno = HashableArray(geno)

            # add probabilities to probability dict
            total_prob += prob
            if geno not in prob_dict:
                prob_dict[geno] = [index, prob]
            else:
                prob_dict[geno][1] += prob

            # add genotype->ID edge
            self._patients_graph.add_edge(geno, patient_id, probability=0, geno_num=index)
            # add subclasses alleles
            subclasses = self.add_subclasses(geno, subclasses)

        # set normalized probabilities to the last patient
        for g, (_, probability) in prob_dict.items():
            self._patients_graph.edges[g, last_patient]['probability'] = \
                probability / total_prob

        return subclasses

    def create_patients_graph_and_find_candidates(self, f_patients: str, verbose: bool) -> dict[int, Sequence[int]]:
        subclasses = self.create_patients_graph(f_patients)
        for subclass in tqdm(subclasses, desc="finding matching candidates", disable=not verbose):
            if self._graph.in_nodes(subclass.subclass):
                self.find_candidates(subclass, self._patients_graph.neighbors(subclass))
        return self.patients

    def score_matches(self, mismatch: int, results_df: pd.DataFrame, donors_info: Iterable[str],
                      patient: int, threshold: float, cutof: int,
                      matched: Set[int]) -> Tuple[Set[int], int, pd.DataFrame]:
        """
        Given a number of mismatches and a patient, this function will return a dictionary
        of all matching donors found in the data with the specific number of mismatches,
        sorted by their probability for a match.

        :param mismatch: number of mismatch to search. could be 0, 1, 2, 3.
        :param results_df: A df storing the matching results.
        :param donors_info: a list of fields from Database to add to the matching results.
        :param patient: patient ID.
        :param threshold: Minimal score value for a valid match. default is 0.1.
        :param cutof: Maximum number of matches to return. default is 50.
        :param matched: A set of donors ID that have already matched for this patient.
        :return: a dictionary of all matching donors and their matching properties.
        """
        if len(matched) >= cutof:
            return matched, 0, results_df

        # a loop that set the scores for all the matching candidates.
        patient_scores = {}
        for hla_id in self._patients_graph.neighbors(patient):
            for i in self._patients_graph.get_edge_data(patient, hla_id)['weight']:
                # match_info = (probability of patient's genotype, number of matches to patient's genotype)
                match_info = self._patients_graph.get_edge_data(patient, hla_id)['weight'][i]
                if match_info[1] == 10 - mismatch:
                    prob = match_info[0]

                    # add the probabilities multiplication of the patient and all the donors that has this genotype
                    # to their matching probabilities.
                    for donor in self.find_donor_from_geno(hla_id):
                        donor_prob = self._graph.get_edge_data(node1=hla_id, node2=donor, node1_id=True)
                        if donor in patient_scores:
                            patient_scores[donor][0] += prob * donor_prob
                            if donor_prob > patient_scores[donor][2]:
                                patient_scores[donor][1:] = [hla_id, donor_prob]

                        else:
                            patient_scores[donor] = [prob * donor_prob, hla_id, donor_prob]

        ids_scores = []
        count_matches = 0

        # sort matching according to their probability
        for donor in patient_scores.keys():
            # do not count or match to an already matched donors.
            if donor in matched or patient_scores[donor][0] < threshold:
                continue

            count_matches += 1
            ids_scores.append((donor, patient_scores[donor][0]))

        ids_scores.sort(reverse=True, key=lambda x: x[1])

        add_donors = {col: [] for col in results_df.columns.values.tolist()}

        # write matching donors to results.
        for donor, score in ids_scores:
            if len(matched) >= cutof:
                break
            matched.add(donor)
            self.append_matching_donor(add_donors, donors_info, patient, donor, score * 100, mismatch)
            # mcg = self.print_most_common_genotype(don_id=donor, pat_geno=self.patients[patient])
            # mcd_match = donor_mismatch_format(self.graph.node_value_from_id(patient_scores[donor][1]),
            #                                   self.patients[patient])
            # donor_info = ",".join([str(v) for v in tuple(DONORS_DB[donor].values())[:-1]])

        results_df = pd.concat([results_df, pd.DataFrame(add_donors)], ignore_index=True)
        return matched, count_matches, results_df

    def append_matching_donor(self, add_donors: Dict, donors_info: Iterable[str],
                              patient: int, donor: int, match_prob: float, mm_number: int) -> None:
        """add a donor to the matches dictionary"""

        add_donors["Patient_ID"].append(patient)
        add_donors["Donor_ID"].append(donor)
        allele_prob = self.probability_to_allele(don_id=donor, pat_geno=self.patients[patient])
        add_donors["Match_Probability_A_1"].append(allele_prob[0])
        add_donors["Match_Probability_A_2"].append(allele_prob[1])
        add_donors["Match_Probability_B_1"].append(allele_prob[2])
        add_donors["Match_Probability_B_2"].append(allele_prob[3])
        add_donors["Match_Probability_C_1"].append(allele_prob[4])
        add_donors["Match_Probability_C_2"].append(allele_prob[5])
        add_donors["Match_Probability_DQB1_1"].append(allele_prob[6])
        add_donors["Match_Probability_DQB1_2"].append(allele_prob[7])
        add_donors["Match_Probability_DRB1_1"].append(allele_prob[8])
        add_donors["Match_Probability_DRB1_2"].append(allele_prob[9])
        add_donors["Matching_Probability"].append(match_prob)
        add_donors["Number_Of_Mismatches"].append(mm_number)
        add_donors["Permissive/Non-Permissive"].append("-")  # TODO: add permissiveness algorithm

        # add the other given fields to the results
        for field in donors_info:
            add_donors[field].append(DONORS_DB.loc[donor, field])
