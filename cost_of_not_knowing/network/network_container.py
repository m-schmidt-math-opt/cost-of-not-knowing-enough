# Author: Martin Schmidt
#         martin.schmidt@uni-trier.de
#         Mathias Sirvent
#         sirvent@posteo.de


import logging
from typing import Optional

import networkx as nx

from coding.cost_of_not_knowing.gaslibparser.gaslibparser.gaslib_parser import GasLibParser
from coding.cost_of_not_knowing.network.pipe_container import PipeContainer

logger = logging.getLogger(__name__.split(".")[-1])


class NetworkContainer:
    """
    The network container model. This class wraps the GasLibParser and executes the methods that parse the data and
    setup the network. Finally it holds the NetworkX object.

    Note that the GasLibParser uses the namedtuble Pipe to hold the data of a Pipe.  As we do not want to edit the
    general GasLibParser, we extend the namedtuple by a PipeContainer to add additional functionality and
    calculations here.
    """

    def __init__(self, net_file, scn_file, sampling_steps: Optional[tuple] = None, is_box_algorithm=False):
        """
        Constructor of the network container model.
        """

        # GasLib parser object
        self._gaslib_parser = GasLibParser(net_file, scn_file)

        # Network (obtained by constructing it using the parsed data)
        self.network = None

        # Parse the data and setup the network
        self.sampling_runtime = None
        self._parse_gaslib_data(sampling_steps, is_box_algorithm)
        self._setup_network()

    def handle_new_point(
        self,
        iteration: int,
        pipe: PipeContainer,
        pressure_in,  # bar
        pressure_out,  # bar
        mass_flow,  # kg_per_second,,
        error_value: float,
        absolute_error: bool,
    ) -> tuple:
        """
        Handles a new sampling point and writes it into the successful points or into the sampling points. Returns a
        boolean, if pipe is optimal and the error value.
        """
        sol_triple = (pressure_in, pressure_out, mass_flow)
        val = abs(pipe.get_f(*sol_triple)) if absolute_error else pipe.get_relative_error(*sol_triple)
        log_helper = "F" if absolute_error else "Relative Error"
        if val >= error_value:
            logger.info(f"Pipe {pipe.arc_id} not optimal. {log_helper} is {val}. Add new sampling point.")
            self._add_sampling_point(
                iteration,
                pipe,
                pressure_in,
                pressure_out,
                mass_flow,
            )
            return False, val
        logger.info(f"Pipe {pipe.arc_id} optimal. {log_helper} is {val}.")
        self._add_successful_point(
            iteration,
            pipe,
            pressure_in,
            pressure_out,
            mass_flow,
        )
        return True, val

    @staticmethod
    def _add_successful_point(
        iteration: int,
        pipe: PipeContainer,
        pressure_in: float,  # bar
        pressure_out: float,  # bar
        mass_flow: float,  # kg_per_second,
    ):
        """
        Adds a successful point to the network.
        """
        pipe.add_successful_point(iteration, pressure_in, pressure_out, mass_flow)

    @staticmethod
    def _add_sampling_point(
        iteration: int,
        pipe: PipeContainer,
        pressure_in: float,  # bar
        pressure_out: float,  # bar
        mass_flow: float,  # kg_per_second,
    ):
        """
        Adds a sampling point to the network.
        """
        pipe.add_sampling_point(iteration, pressure_in, pressure_out, mass_flow, True)

    def _parse_gaslib_data(
        self,
        sampling_steps: Optional[tuple],
        is_box_algorithm: bool,
    ):
        """
        Parsing the given GasLib data.
        """

        # Parsing
        self._gaslib_parser.parse()

        # Extracting node data
        self._entries = self._gaslib_parser.entry_data
        self._exits = self._gaslib_parser.exit_data
        self._innodes = self._gaslib_parser.innode_data
        self._nodes = {**self._entries, **self._exits, **self._innodes}

        # Extracting arc data
        # Use pipe container as extension of the Pipe namedtuple to store data and implement additional methods
        _pipes = {}
        sampling_runtime = 0.0
        for pipe_id, pipe in self._gaslib_parser.pipe_data.items():
            pipe = PipeContainer(
                pipe,
                self._nodes.get(pipe.from_node),
                self._nodes.get(pipe.to_node),
                sampling_steps,
                is_box_algorithm,
            )
            sampling_runtime += pipe.sampling_runtime
            _pipes.update({pipe_id: pipe})
        self.sampling_runtime = sampling_runtime if sampling_steps else 0.0
        self._pipes = _pipes
        self._compressor_stations = self._gaslib_parser.compressor_station_data

        logger.info("Parsing ... done")

    def _setup_network(self):
        """
        Constructing the NetworkX graph.
        """
        self.network = nx.MultiDiGraph(name=self._gaslib_parser.net_name)

        # Add nodes
        for nodeId, nodeData in self._entries.items():
            self.network.add_node(nodeId, data=nodeData)
        for nodeId, nodeData in self._exits.items():
            self.network.add_node(nodeId, data=nodeData)
        for nodeId, nodeData in self._innodes.items():
            self.network.add_node(nodeId, data=nodeData)

        # Add arcs
        for arcId, arcData in self._pipes.items():
            self.network.add_edge(arcData.from_node, arcData.to_node, data=arcData)
        for arcId, arcData in self._compressor_stations.items():
            self.network.add_edge(arcData.from_node, arcData.to_node, data=arcData)
        logger.info("Constructing the network ... done")
