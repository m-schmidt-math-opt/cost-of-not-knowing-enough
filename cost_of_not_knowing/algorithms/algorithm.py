from __future__ import annotations

import json
import time
from abc import abstractmethod
from pathlib import Path

from coding.cost_of_not_knowing.gaslibparser.gaslibparser.unit_helper import unit


class Algorithm:
    """
    The algorithm class.
    """

    def __init__(self, index: int, net_file: str):
        """
        Constructor of the algorithm model class.
        """
        testset = Path(net_file).parent.stem
        timestr = time.strftime("%Y_%m_%d_%H_%M_%S")
        self._result_directory = Path(f"results/{self.__class__.__name__}/{timestr}_{index}_{testset}/")
        self._result_directory.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def execute(self):
        """
        Entry point to execute the algorithm.
        """

        pass

    def _write_solution(
        self,
        node_to_pressure: dict,
        arcs_to_flows: dict,
        bnd_nodes_to_flows: dict,
        compressors_to_press_inc: dict,
    ):
        """
        Helper to write out solutions the same way for all algorithms.
        """
        self._write_single_solution("node_to_pressure", node_to_pressure)
        self._write_single_solution("arcs_to_flows", arcs_to_flows)
        self._write_single_solution("bnd_nodes_to_flows", bnd_nodes_to_flows)
        self._write_single_solution("compressors_to_press_inc", compressors_to_press_inc)

    def _write_single_solution(
        self,
        name: str,
        dic: dict,
    ):
        """
        Helper to write out solution.
        """
        new_dict = {}
        for key, value in dic.items():
            if isinstance(value, unit.Quantity):
                value = value.to_tuple()
            if isinstance(key, tuple):
                key = "_".join(key)
            new_dict.update({key: value})
        with open(str(self._result_directory.joinpath(f"{name}.json")), "w") as f:
            f.write(json.dumps(new_dict, indent=4))
