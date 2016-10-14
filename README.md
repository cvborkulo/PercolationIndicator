PercolationIndicator
==========

This method captures the dynamics of a network in a single parameter: the percolation indicator. The indicator is computed by combining intensive longitudinal binary data of a system with its network structure (with temporal relationships). The dynamics of the network is modeled by two independent Poisson processes. The ratio of the parameters of these processes (infection rate and recovery rate) is used to obtain an indication of the development of the network over time. When the percolation indicator pi <= 1, activity in the network will die out. When pi > 1, the network will remain infected indefinitely. Can deal with binary data and an unweighted adjacency matrix!
