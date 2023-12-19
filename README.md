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
The damage probabilities for each type of component are specified as follows. 

Note that the partial damage (i.e., part

# Outputs
Both scripts produce two outputs that define the network:
* ./outputs/nodes/nodes.shp: nodes shapefile that includes the following attributes
* ./outputs/edges/edges.shp: edges shapefile that includes the following attributes

# Reference
 Dunton, A. (2023). Model selection for risk analysis of wastewater networks. https://doi.org/10.48550/arXiv.2312.06623
