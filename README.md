# ww-model-selection
This repository includes code for simulating the effect of component damage on the functionality of wastewater networks. Specifically, two methods are used for assessing the functionality, with two corresponding scripts here:
* flow.py: uses flow analysis (specifically SWMM) for functionality assessment
* connectivity.py: uses connectivity analysis for functionality assessment

For each of these cases, a performance measure (PM) is used to describe the loss of functionality to each node. PM ranges from 0 to 1, where 1 is no loss of functionality and 0 is complete loss of functionality. See Dunton (2023) for the formulation of PM for connectivity analysis and flow analysis. 

Each node has a deterministic value of PM for each simulation of the damage. That is, the simulations are mapping the uncertainty in the damage state of components to the uncertainty in the value of the nodal PM. To determine the number of simulations to use, a convergence criterion is used. Here, the simulations continue until the average coefficient of variation (COV) of the estimate of PM over the nodes is less than 0.1. See Dunton (2023) for more information. The convergence criteria is checked every 100 simulations, and some information is printed to the terminal. Once the convergence criteria is met, the total simulation time and final average COV of PM are also printed to the terminal.

# Specifying the Input Model
The analysis starts from an existing SWMM model representing the undamaged network. This model should be included in ./inputs/, and the path to this model must be specified in connectivity.py or flow.py as INP_PATH. This repository includes 4 example SWMM models that were used to generate the results in Dunton (2023). By running simulations for 2 model granularities and 3 model fidelities, Dunton (2023) assesses the necessary granularity and fidelity for accurate risk analysis of wastewater networks. The following table summarizes the experimental design and how the files in this repository were used to generate the results:

**Summary of Experimental Design from Dunton (2023)**
| Model Name | Granularity | Fidelity | Script Used | Undamaged Model |
|---|---|---|---|---|
|DF|Full Network|Flow Analysis, Dynamic Wave Routing|flow.py|./inputs/full_dynwave|
|KF|Full Network|Flow Analysis, Kinematic Wave Routing|flow.py|./inputs/full_kinwave|
|CF|Full Network|Connectivity Analysis|connectivity.py|./inputs/full_kinwave **or** ./inputs/full_dynwave|
|DA|Arterial Network|Flow Analysis, Dynamic Wave Routing|flow.py|./inputs/arterial_dynwave|
|KA|Arterial Network|Flow Analysis, Kinematic Wave Routing|flow.py|./inputs/arterial_kinwave|
|CA|Arterial Network|Connectivity Analysis|connectivity.py|./inputs/arterial_kinwave **or** ./inputs/arterial_dynwave|

# Specifying the Damage Probabilities
The damage probabilities for each type of component are specified as inputs in connectivity.py and flow.py. For connectivity.py, components are either considered fully functionally or fully damaged. As such, the failure probabilities of pumps, tanks, and pipes are specified as follows:
* PF_PUMP = 0.1
* PF_TANK = 0.1
* PF_PIPE = 0.01

Flow analysis allows us to additionally model partial damage to the components. Specifically, we model partial damage to the pipes by reducing their flow capacity (i.e., diameter). In the flow script, the probability of damage to the pipes is represented by PF_PIPE. PF_PUMP and PF_TANK are the same as above. For flow.py, PF_PIPE is a list of tuples, where each tuple in the list represents a damage state. The first value in each tuple is the fraction of the flow that is able to pass through the pipe in the corresponding damage state; the second value in each tuple is the probability of the corresponding damage state.

# Outputs
Both connectivity.py and flow.py produce two outputs that define the network:
* ./outputs/nodes/nodes.shp: nodes shapefile that includes the following attributes for each node
    * node: node index
    * N: number of simulations completed
    * pm_avg: average performance measure over the simulations
    * pm_200, pm_201, ..., pm_299: performance measures for the last 100 simulations 
* ./outputs/edges/edges.shp: edges shapefile that includes the following attributes for each edge
    * from: node index of the source node of the edge
    * to: node index of the destination node of the edge

# Reference
 Dunton, A. (2023). Model selection for risk analysis of wastewater networks. https://doi.org/10.48550/arXiv.2312.06623
