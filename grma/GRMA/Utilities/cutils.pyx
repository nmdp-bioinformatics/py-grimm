#cython: language_level=3

import cython
import numpy as np
cimport numpy as np

ctypedef np.uint8_t UINT8
ctypedef np.int8_t INT8
ctypedef np.uint16_t UINT16
ctypedef np.uint32_t UINT32



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[UINT32, ndim=2] cdrop_less_than_7_matches(np.ndarray[UINT32, ndim=1] ids,
                                                         np.ndarray[INT8, ndim=1] similarities):
    cdef:
        np.ndarray[UINT32, ndim=2] ret
        UINT32 count, i

    ret = np.zeros((len(ids), 2), dtype=np.uint32)
    count = 0
    for i in range(len(ids)):
        if similarities[i] != -1:
            ret[count, 0] = ids[i]
            ret[count, 1] = similarities[i]
            count += 1
    return ret[:count, :]


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[INT8, ndim=1] ccheck_similarity(np.ndarray[UINT16, ndim=1] patients_geno,
                                                 np.ndarray[UINT16, ndim=2] donors_genos,
                                                 np.ndarray[UINT8, ndim=1] allele_range, UINT8 init_count_similar):
    """
    Takes 2 genotypes to check similarity between \n
    A list of even numbers between 0-10 of alleles to check. \n
    And the number of alleles that are known to match.
    Count allele matches between the genotypes.

    :returns: the number of similarities between the given genotypes.
    If the count is bigger than 7, it is returned, else returns -1.
    """
    cdef:
        np.ndarray[INT8, ndim=1] similarities
        np.ndarray[UINT16, ndim=1] donors_geno
        UINT8 count_similar, counted, number_of_alleles, allele_num
        UINT16 patient_alleles0, patient_alleles1, donor_alleles0, donor_alleles1

    similarities = np.zeros(len(donors_genos), dtype=np.int8)
    number_of_alleles = len(allele_range)
    for i in range(len(donors_genos)):
        counted: int = init_count_similar
        count_similar: int = init_count_similar
        donors_geno = donors_genos[i]
        for j in range(number_of_alleles):
            allele_num = allele_range[j]
            counted += 2
            patient_alleles0 = patients_geno[allele_num]
            patient_alleles1 = patients_geno[allele_num + 1]
            donor_alleles0 = donors_geno[allele_num]
            donor_alleles1 = donors_geno[allele_num + 1]

            if patient_alleles0 == donor_alleles0:
                if patient_alleles1 == donor_alleles1:
                    count_similar += 2
                else:
                    count_similar += 1

            elif patient_alleles1 == donor_alleles1:
                count_similar += 1

            elif patient_alleles0 == donor_alleles1:
                count_similar += 1

            elif patient_alleles1 == donor_alleles0:
                count_similar += 1

            if counted - count_similar > 3:
                similarities[i] = -1
                break
        if 10 - count_similar > 3:
            similarities[i] = -1
        else:
            similarities[i] = count_similar

    return similarities


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef UINT32 chash(np.ndarray[UINT16, ndim=1] arr):
    cdef UINT32 h = 17
    cdef UINT8 i
    for i in range(len(arr)):
        h = h * 31 + arr[i]
    return h
