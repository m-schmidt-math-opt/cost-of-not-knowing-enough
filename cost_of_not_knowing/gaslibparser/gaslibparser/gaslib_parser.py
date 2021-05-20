import re
from xml.etree import cElementTree

import numpy as np

from .connections import Pipe, CompressorStation, Valve, ControlValve, ShortPipe, Resistor
from .nodes import Node, Entry, Exit
from .unit_helper import unit


class GasLibParser(object):
    """Class for parsing GasLib tests."""

    def __init__(self, gaslib_net_file, gaslib_scn_file):
        """Constructor."""
        self.gaslib_net_file = gaslib_net_file
        self.gaslib_scn_file = gaslib_scn_file
        self.namespaces = {"framework": "http://gaslib.zib.de/Framework", "gas": "http://gaslib.zib.de/Gas"}
        self.net_name = None
        # Nodes
        self.entry_data = {}
        self.exit_data = {}
        self.innode_data = {}
        # Arcs ; todo really all needed
        self.pipe_data = {}
        self.compressor_station_data = {}
        self.valve_data = {}
        self.control_valve_data = {}
        self.short_pipe_data = {}
        self.resistor_data = {}

    def parse(self):
        """Main parsing method."""
        self._parse_net_file()
        self._parse_scn_file()

    def _parse_net_file(self):
        etree = cElementTree.parse(self.gaslib_net_file)
        self.net_name = etree.find(".//framework:title", self.namespaces).text
        nodes_element = etree.find("framework:nodes", self.namespaces)
        self._parse_nodes(nodes_element)
        connections_element = etree.find("framework:connections", self.namespaces)
        self._parse_connections(connections_element)

    def _parse_nodes(self, nodes_element):
        """parse the nodes section"""
        sources_iter = nodes_element.iterfind("gas:source", self.namespaces)
        self._parse_sources(sources_iter)
        sinks_iter = nodes_element.iterfind("gas:sink", self.namespaces)
        self._parse_sinks(sinks_iter)
        innodes_iter = nodes_element.iterfind("gas:innode", self.namespaces)
        self._parse_innodes(innodes_iter)

    def _parse_innodes(self, innodes_iter):
        for innode in innodes_iter:
            node_dict = {"node_id": innode.get("id"), "pos": [float(innode.get("x")), float(innode.get("y"))]}
            for element in innode.iter():
                tag = self.strip_namespace(element)
                # skip the innode tags
                if tag == "innode":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                node_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            i_data = Node(**node_dict)
            self.innode_data[i_data.node_id] = i_data

    def _parse_sinks(self, sinks_iter):
        for sink in sinks_iter:
            exit_dict = {"node_id": sink.get("id"), "pos": [float(sink.get("x")), float(sink.get("y"))]}
            for element in sink.iter():
                tag = self.strip_namespace(element)
                if tag == "sink":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                exit_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            e_data = Exit(**exit_dict)
            self.exit_data[e_data.node_id] = e_data

    def _parse_sources(self, sources_iter):
        for source in sources_iter:
            entry_dict = {"node_id": source.get("id"), "pos": [float(source.get("x")), float(source.get("y"))]}
            for element in source.iter():
                tag = self.strip_namespace(element)
                if tag == "coefficient-A-heatCapacity":
                    new_name = "heat_coeff_A"
                elif tag == "coefficient-B-heatCapacity":
                    new_name = "heat_coeff_B"
                elif tag == "coefficient-C-heatCapacity":
                    new_name = "heat_coeff_C"
                elif tag == "source":
                    continue
                elif tag == "gasTemperature":
                    new_name = "gas_temp"
                else:
                    new_name = self.cc_to_us(tag)
                entry_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            e_data = Entry(**entry_dict)
            self.entry_data[e_data.node_id] = e_data

    def _parse_connections(self, connections_element):
        """parse the connections section"""
        pipes_iter = connections_element.iterfind("gas:pipe", self.namespaces)
        self._parse_pipes(pipes_iter)
        compressors_iter = connections_element.iterfind("gas:compressorStation", self.namespaces)
        self._parse_compressor_stations(compressors_iter)
        cvalves_iter = connections_element.iterfind("gas:controlValve", self.namespaces)
        self._parse_control_valves(cvalves_iter)
        valves_iter = connections_element.iterfind("gas:valve", self.namespaces)
        self._parse_valves(valves_iter)
        spipes_iter = connections_element.iterfind("gas:shortPipe", self.namespaces)
        self._parse_short_pipes(spipes_iter)
        resistors_iter = connections_element.iterfind("gas:resistor", self.namespaces)
        self._parse_resistors(resistors_iter)

    def _parse_compressor_stations(self, compressors_iter):
        for compressor in compressors_iter:
            compressor_dict = {
                "arc_id": compressor.get("id"),
                "from_node": compressor.get("from"),
                "to_node": compressor.get("to"),
            }  # self.init_compressor()
            for element in compressor.iter():
                tag = self.strip_namespace(element)
                if tag == "compressorStation":
                    continue
                elif tag == "dragFactorIn" or tag == "dragFactorOut":
                    continue
                elif tag == "diameterIn" or tag == "diameterOut":
                    continue
                elif tag == "pressureLossIn" or tag == "pressureLossOut":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                compressor_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            cs_data = CompressorStation(**compressor_dict)
            self.compressor_station_data[cs_data.arc_id] = cs_data

    def _parse_control_valves(self, cvalves_iter):
        for cvalve in cvalves_iter:
            cvalve_dict = {"arc_id": cvalve.get("id"), "from_node": cvalve.get("from"), "to_node": cvalve.get("to")}
            for element in cvalve.iter():
                tag = self.strip_namespace(element)
                if tag == "controlValve":
                    continue
                elif tag == "dragFactorIn" or tag == "dragFactorOut":
                    continue
                elif tag == "diameterIn" or tag == "diameterOut":
                    continue
                elif tag == "pressureLossIn" or tag == "pressureLossOut":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                cvalve_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            cv_data = ControlValve(**cvalve_dict)
            self.control_valve_data[cv_data.arc_id] = cv_data

    def _parse_valves(self, valves_iter):
        for valve in valves_iter:
            valve_dict = {"arc_id": valve.get("id"), "from_node": valve.get("from"), "to_node": valve.get("to")}
            for element in valve.iter():
                tag = self.strip_namespace(element)
                if tag == "valve":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                valve_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            valve_data = Valve(**valve_dict)
            self.valve_data[valve_data.arc_id] = valve_data

    def _parse_pipes(self, pipes_iter):
        for pipe in pipes_iter:
            pipe_dict = {"arc_id": pipe.get("id"), "from_node": pipe.get("from"), "to_node": pipe.get("to")}
            for element in pipe.iter():
                tag = self.strip_namespace(element)
                if tag == "pipe":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                pipe_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            p_data = Pipe(**pipe_dict)
            self.pipe_data[p_data.arc_id] = p_data

    def _parse_short_pipes(self, spipes_iter):
        for spipe in spipes_iter:
            spipe_dict = {"arc_id": spipe.get("id"), "from_node": spipe.get("from"), "to_node": spipe.get("to")}
            for element in spipe.iter():
                tag = self.strip_namespace(element)
                if tag == "shortPipe":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                spipe_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            sp_data = ShortPipe(**spipe_dict)
            self.pipe_data[sp_data.arc_id] = sp_data

    def _parse_resistors(self, resistors_iter):
        for resistor in resistors_iter:
            resistor_dict = {
                "arc_id": resistor.get("id"),
                "from_node": resistor.get("from"),
                "to_node": resistor.get("to"),
            }
            for element in resistor.iter():
                tag = self.strip_namespace(element)
                if tag == "resistor":
                    continue
                else:
                    new_name = self.cc_to_us(tag)
                resistor_dict[new_name] = self._parse_value(element.get("value"), element.get("unit"))
            r_data = Resistor(**resistor_dict)
            self.resistor_data[r_data.arc_id] = r_data

    def _parse_scn_file(self):
        etree = cElementTree.parse(self.gaslib_scn_file)
        nodes_iter = etree.iterfind(".//gas:node", self.namespaces)
        for node_element in nodes_iter:
            self._parse_scn_node(node_element)

    def _parse_scn_node(self, node_element):
        scn_intervals = {
            "pressure": (unit.Quantity(0, "bar"), unit.Quantity(np.Infinity, "bar")),
            "flow": (unit.Quantity(0, "m_cube_per_hour"), unit.Quantity(np.Infinity, "m_cube_per_hour")),
        }
        node_id = node_element.get("id")
        node_type = node_element.get("type")
        for bound in node_element.iter():
            # skip the node xml lines
            if self.strip_namespace(bound) == "node":
                continue
            # pressure or flow
            bound_type = self.strip_namespace(bound)
            # lower, upper or both
            bound_sense = bound.get("bound")
            # create bound value with the corresponding unit
            bound_value = self._parse_value(bound.get("value"), bound.get("unit"))
            # update bounds arrconding to bound_sense
            if bound_sense == "lower" or bound_sense == "both":
                s_max = scn_intervals[bound_type][1]
                scn_intervals[bound_type] = (bound_value, s_max)
            if bound_sense == "upper" or bound_sense == "both":
                s_min = scn_intervals[bound_type][0]
                scn_intervals[bound_type] = (s_min, bound_value)

        if node_type == "exit":
            node_data = self.exit_data[node_id]._asdict()
        elif node_type == "entry":
            node_data = self.entry_data[node_id]._asdict()
        else:
            node_data = None

        node_data["nomination_min"] = scn_intervals["flow"][0]
        node_data["nomination_max"] = scn_intervals["flow"][1]

        net_pressure_interval = (node_data["pressure_min"], node_data["pressure_max"])
        result_pressure_interval = self._intersect_intervals(scn_intervals["pressure"], net_pressure_interval)

        # convert pressures to bar to avoid using bar gauge units
        node_data["pressure_min"] = result_pressure_interval[0].to("bar")
        node_data["pressure_max"] = result_pressure_interval[1].to("bar")

        if node_type == "exit":
            self.exit_data[node_id] = Exit(**node_data)
        elif node_type == "entry":
            self.entry_data[node_id] = Entry(**node_data)
        else:
            # there is additional info on inner nodes
            pass

    @staticmethod
    def _intersect_intervals(a_int, b_int):
        (amin, amax) = a_int
        (bmin, bmax) = b_int
        return max(amin, bmin), min(amax, bmax)

    @staticmethod
    def _parse_value(value_string, unit_string):
        val = np.float64(value_string)
        if unit_string is not None:
            unit_numeric = re.sub(r"\D.*", "", unit_string)
            factor = np.float64(1 if unit_numeric == "" else unit_numeric)
            unit_alpha = re.findall(r"\D.*", unit_string)[0]

            # pint does not handle offset units very well, thus do it here
            if unit_alpha == "barg":
                unit_alpha = "bar"
                val += 1.01325
        else:
            factor = 1
            unit_alpha = None
        return unit.Quantity(val * factor, unit_alpha)

    @staticmethod
    def strip_namespace(element):
        return re.sub("[^}]*}(.*)", r"\1", element.tag)

    @staticmethod
    def cc_to_us(name):
        """"""
        # taken from
        # http://stackoverflow.com/questions/1175208/elegant-python-function-to-convert-camelcase-to-camel-case
        s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
        return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()

    # # Where does this come from?
    # def init_compressor(self):
    #     """"""
    #     if self.net_name == "GasLib_40":
    #         compressor_dict = {"max_pressure_increase": 24.33131262,
    #                            "max_pressure_ratio": 1.521214713,
    #                            "min_pressure_increase": 1.883212208,
    #                            "min_pressure_ratio": 1.060722826}
    #     elif self.net_name == "GasLib_135":
    #         compressor_dict = {"max_pressure_increase": 31.70737918,
    #                            "max_pressure_ratio": 1.806683035,
    #                            "min_pressure_increase": 1.034326951,
    #                            "min_pressure_ratio": 1.033351131}
    #     # using a generic compressor dict for now
    #     else:
    #         compressor_dict = {"max_pressure_increase": 31.70737918,
    #                            "max_pressure_ratio": 1.806683035,
    #                            "min_pressure_increase": 1.034326951,
    #                            "min_pressure_ratio": 1.033351131}
    #     return compressor_dict
