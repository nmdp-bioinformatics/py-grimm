from __future__ import annotations

import pickle
from os import PathLike
from typing import Union

import numpy as np

from grma.utilities.geno_representation import HashableArray
from grma.match.lol_graph import LolGraph

NODES_TYPES = Union[int, HashableArray]


class Graph(object):
    """Graph wrapper class for LOLGraph"""
    __slots__ = "_map_node_to_number", "_graph"

    def __init__(self, lol_properties: dict):
        self._map_node_to_number = lol_properties["map_node_to_number"]

        self._graph = LolGraph(index_list=lol_properties["index_list"],
                               neighbors_list=lol_properties["neighbors_list"],
                               weights_list=lol_properties["weights_list"],
                               map_number_to_num_node=lol_properties["map_number_to_num_node"],
                               map_number_to_arr_node=lol_properties["map_number_to_arr_node"],
                               arrays_start=lol_properties["arrays_start"],
                               directed=lol_properties["directed"],
                               weighted=lol_properties["weighted"])

    def in_nodes(self, node: NODES_TYPES) -> bool:
        """return True if the given node is in the graph and false otherwise"""
        return node in self._map_node_to_number

    def get_node_id(self, node: NODES_TYPES) -> bool:
        """return True if the given node is in the graph and false otherwise"""
        return self._map_node_to_number.get(node, None)

    def get_edge_data(self, node1: NODES_TYPES, node2: NODES_TYPES,
                      node1_id: bool = False, node2_id: bool = False, default: object = None):
        """
        Returns edge properties. If does not exist return default
        """
        exception_val = -1
        node1_num = self._map_node_to_number[node1] if not node1_id else node1
        node2_num = self._map_node_to_number[node2] if not node2_id else node2
        ret = self._graph.get_edge_data(node1_num, node2_num)
        return default if ret == exception_val else ret

    def class_neighbors(self, node: NODES_TYPES | int, search_lol_id: bool = False):
        node_num = self._map_node_to_number[node] if not search_lol_id else node
        neighbors_list = self._graph.neighbors_unweighted(node_num)

        neighbors_list_values = np.ndarray([len(neighbors_list), 10], dtype=np.uint16)
        for i, neighbor in enumerate(neighbors_list):
            neighbors_list_values[i, :] = self._graph.arr_node_value_from_id(neighbor)

        return neighbors_list, neighbors_list_values

    def neighbors_unweighted(self, node: NODES_TYPES | int, search_lol_id: bool = False):
        """
        Get node's neighbors. There are some options for what to get (see parameters).
         - Search node by its name: search_lol_id=False
         - Search node by its lol ID: search_lol_id=True
        :return: tuple of lists: (neighbor's IDs, neighbor's values, weights)
        """
        node_num = self._map_node_to_number[node] if not search_lol_id else node
        neighbors_list = self._graph.neighbors_unweighted(node_num)

        neighbors_list_values = [0] * len(neighbors_list)
        for i, neighbor in enumerate(neighbors_list):
            neighbors_list_values[i] = self._graph.arr_node_value_from_id(neighbor)

        return neighbors_list, neighbors_list_values

    def neighbors(self, node: NODES_TYPES | int, search_lol_id: bool = False) -> zip[tuple[int, float]] | list[int]:
        """
        Get node's neighbors. There are some options for what to get (see parameters).
         - Search node by its name: search_lol_id=False
         - Search node by its lol ID: search_lol_id=True
        :return: tuple of lists: (neighbor's IDs, neighbor's values, weights)
        """
        node_num = self._map_node_to_number[node] if not search_lol_id else node
        if self._graph.is_weighted():
            neighbors_list, weights_list = self._graph.neighbors_weighted(node_num)
        else:
            neighbors_list = self._graph.neighbors_unweighted(node_num)

        neighbors_list_values = [0] * len(neighbors_list)
        for i, neighbor in enumerate(neighbors_list):
            neighbors_list_values[i] = self.node_value_from_id(neighbor)

        if self._graph.is_weighted():
            return zip(*(neighbors_list_values, weights_list))
        return neighbors_list_values

    def neighbors_2nd(self, node):
        node_num = self._map_node_to_number[node]
        r0, r1 = self._graph.neighbors_2nd(node_num)
        return r0[:-1], r1

    def node_value_from_id(self, node_id: int) -> NODES_TYPES:
        """convert lol ID to node value"""
        if node_id < self._graph.array_start:
            return self._graph.num_node_value_from_id(node_id)
        return self._graph.arr_node_value_from_id(node_id)

    @classmethod
    def from_pickle(cls, path: Union[str, PathLike]):
        graph_dict = pickle.load(open(path, "rb"))
        return cls(graph_dict)
