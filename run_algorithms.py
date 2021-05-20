import logging

from coding.cost_of_not_knowing.algorithms.benchmark import Benchmark, Solver
from coding.cost_of_not_knowing.algorithms.box_algorithm import BoxAlgorithm

logging.basicConfig(level=logging.INFO)

trees = [
    "GasLib-4-Tree",
    "GasLib-11",
]
circles = [
    "GasLib-4",
    "GasLib-24",
]

# Trees
for test_set in trees:
    net_file = f"tests/data/{test_set}/net.net"
    scn_file = f"tests/data/{test_set}/scn.scn"
    algorithms = Benchmark(0, net_file, scn_file, Solver.scip)
    algorithms.execute()
    conf = [None, (10.0, 1.0), (5.0, 1.0), (2.5, 1.0), (1.0, 1.0), (0.5, 1.0), (0.1, 1.0)]
    for i, sp in enumerate(conf):
        box_algorithm = BoxAlgorithm(i, net_file, scn_file, sampling_steps=sp, error_value=0.01, absolute_error=False)
        box_algorithm.execute()

# Circles
for test_set in circles:
    net_file = f"tests/data/{test_set}/net.net"
    scn_file = f"tests/data/{test_set}/scn.scn"
    algorithms = Benchmark(0, net_file, scn_file, Solver.scip)
    algorithms.execute()
    conf = [None, (10.0, 100.0), (5.0, 50.0), (5.0, 10.0), (5.0, 5.0), (1.0, 50.0), (1.0, 10.0), (1.0, 5.0), (1.0, 1.0)]
    for i, sp in enumerate(conf):
        box_algorithm = BoxAlgorithm(i, net_file, scn_file, sampling_steps=sp, error_value=0.01, absolute_error=False)
        box_algorithm.execute()
