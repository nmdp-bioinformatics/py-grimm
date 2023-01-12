from typing import List, Dict, Set

import numpy as np
from tqdm import tqdm
import gc
from collections import OrderedDict

from grma.donorsgraph import Edge
from grma.utilities.utils import print_time, tuple_geno_to_int


class LolBuilder:
    __slots__ = '_directed', '_weighted', '_verbose', '_properties', '_graph'

    def __init__(self, directed: bool, weighted: bool, verbose: bool = False):
        self._directed = directed
        self._weighted = weighted
        self._verbose = verbose
        self._properties = {
            "weighted": self._weighted,
            "directed": self._directed
        }
        self._graph: List[Edge] = []

    def build(self, edge_list: List[Edge], layers: Dict[str, Set]):
        self._graph: List[Edge] = edge_list
        subclasses_start = self._convert(layers)

        neighbors_list, weights_list = self._sort_all(index_list=self._properties["index_list"],
                                                      neighbors_list=self._properties["neighbors_list"],
                                                      weights_list=self._properties["weights_list"])
        if self._weighted:
            weights_list = self._dist_2nd_weights(subclasses_start, self._properties["arrays_start"],
                                                  index_list=self._properties["index_list"],
                                                  neighbors_list=neighbors_list,
                                                  weights_list=weights_list)

        self._properties["neighbors_list"] = neighbors_list
        self._properties["weights_list"] = weights_list

        print_time("Finished creating the lol-matching graph")
        return self._properties

    def _convert(self, layers: Dict[str, Set]):
        """
        convert edgelist to LOL dict-representation according to the donors' graph architecture.
        :param layers: a dictionary of the graph's layers. each layer is a set of all the nodes it contains.
        """
        free = 0
        print_time('(1/6) Convert edgelist to maps from nodes to internal numbers')

        # maps nodes' original value to its lol id.
        # maps only ids and subclasses
        map_node_to_number = OrderedDict()

        for idd in tqdm(layers["ID"], desc="(1.1) Map nodes to internal numbers", disable=not self._verbose):
            map_node_to_number[idd] = free
            free += 1

        subclasses_start = free  # a flag for where the subclasses mapping starts.
        for sub in tqdm(layers["SUBCLASS"], desc="(1.2) Map nodes to internal numbers", disable=not self._verbose):
            map_node_to_number[sub] = free
            free += 1

        # a flag for where the arrays mapping starts.
        arrays_start = free
        # map lol-ids to arrays
        # given an lol_id, the mapping will be map_number_to_arr_node[lol_id - arrays_start, :]
        map_number_to_arr_node = np.zeros((len(layers["GENOTYPE"]), 10), dtype=np.uint16)
        for i, geno in tqdm(enumerate(layers["GENOTYPE"]), desc="(1.3) Map nodes to internal numbers",
                            disable=not self._verbose):
            map_node_to_number[geno] = free
            map_number_to_arr_node[i, :] = geno.np()
            free += 1

        # map classes to lol-id.
        for clss in tqdm(layers["CLASS"], desc="(1.4) Map nodes to internal numbers", disable=not self._verbose):
            # Here also add dictionary {class: id_in_graph}
            map_node_to_number[clss] = free
            free += 1

        print_time('(2/6) Create the opposite map')
        # This map is for the donors Ids only.
        map_number_to_num_node = np.array(
            [x for x, y in tqdm(map_node_to_number.items(), desc="(2/6) Create the opposite map",
                                disable=not self._verbose) if y < subclasses_start], dtype=np.uint32)

        print_time('(3/6) Create the index list')
        num_of_neighbors = np.zeros(len(map_node_to_number))  # set the number of neighbors for every node in

        for edge in tqdm(self._graph, desc="(3.1) Create the index list", disable=not self._verbose):
            num_of_neighbors[map_node_to_number[edge.node1]] += 1

            if not self._directed:
                if edge.node1 != edge.node2:
                    num_of_neighbors[map_node_to_number[edge.node2]] += 1

        index_list = np.zeros(len(num_of_neighbors) + 1, dtype=np.uint32)
        for j in tqdm(range(1, len(num_of_neighbors) + 1), desc="(3.2) Create the index list",
                      disable=not self._verbose):
            index_list[j] = index_list[j - 1] + num_of_neighbors[j - 1]

        del num_of_neighbors
        gc.collect()

        print_time('(4/6) Create the neighbors list')
        if self._directed:
            neighbors_list = -np.ones(len(self._graph))
            if self._weighted:
                weights_list = np.zeros(len(self._graph))
        else:
            neighbors_list = -np.ones(2 * len(self._graph))
            if self._weighted:
                weights_list = np.zeros(2 * len(self._graph))

        space = -np.ones(free, dtype=np.int32)

        for edge in tqdm(self._graph, desc="(4/6) Create the neighbors list", disable=not self._verbose):
            left = map_node_to_number[edge.node1]
            right = map_node_to_number[edge.node2]
            if self._weighted:
                try:
                    weight = float(edge.weight)
                except [IndexError, TypeError]:
                    weight = None

            self._add_weights_and_neighbors(left, right, space, neighbors_list, index_list,
                                            weight=weight, weights_list=weights_list)

            if not self._directed and left != right:
                self._add_weights_and_neighbors(right, left, space, neighbors_list, index_list,
                                                weight=weight, weights_list=weights_list)

        # replace geno hashable array to more efficient representation.
        for array_geno in layers["GENOTYPE"]:
            int_geno = tuple_geno_to_int(array_geno)
            map_node_to_number[int_geno] = map_node_to_number[array_geno]
            del map_node_to_number[array_geno]

        del self._graph
        del layers
        gc.collect()

        # set the lol-properties dictionary
        self._properties["index_list"] = index_list
        self._properties["neighbors_list"] = neighbors_list
        self._properties["weights_list"] = weights_list
        self._properties["map_node_to_number"] = map_node_to_number
        self._properties["map_number_to_num_node"] = map_number_to_num_node
        self._properties["map_number_to_arr_node"] = map_number_to_arr_node
        self._properties["arrays_start"] = arrays_start

        return subclasses_start

    def _add_weights_and_neighbors(self, node1, node2, space, neighbors_list, index_list,
                                   *, weight=None, weights_list=None):
        if space[node1] != -1:
            space[node1] += 1
            i = space[node1]
        else:
            i = index_list[node1]
            space[node1] = i
        neighbors_list[i] = node2
        if self._weighted:
            weights_list[i] = weight

    def _dist_2nd_weights(self, subs_start, subs_end, index_list, neighbors_list, weights_list):
        """
        For each subclass node, add to its weights the number of genotypes connected to it.
        this information is written only in the first node in order that is connected to the subclass.
        """
        print_time("(6/6) Add weights")
        for i in tqdm(range(subs_start, subs_end), desc="(6/6) Add weights", disable=not self._verbose):
            neigh_2nd_total = 0
            start = index_list[i]
            end = index_list[i + 1]

            for n in neighbors_list[start: end]:
                start_n = index_list[n]
                end_n = index_list[n + 1]
                neigh_2nd_total += end_n - start_n

            weights_list[start] = neigh_2nd_total
        return weights_list

    def _sort_all(self, index_list=None, neighbors_list=None, weights_list=None):
        """sort the neighbors for each node"""
        print_time("(5/6) Sort")
        neighbors_list = np.array(neighbors_list, dtype=np.uint32)
        weights_list = np.array(weights_list, dtype=np.float32)

        for number in tqdm(range(len(index_list) - 1), desc="(5/6) Sort", disable=not self._verbose):
            start = index_list[number]
            end = index_list[number + 1]
            neighbors_list[start: end], weights_list[start: end] = self._sort_neighbors(neighbors_list[start: end],
                                                                                        weights_list[start: end])
        return neighbors_list, weights_list

    def _sort_neighbors(self, neighbors_list=None, weights_list=None):
        if self._weighted:
            # sort both neighbors_list and weights_list according to the neighbors_list
            index = neighbors_list.argsort()
            return neighbors_list[index], weights_list[index]
        else:
            return sorted(neighbors_list), weights_list
