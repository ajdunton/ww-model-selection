# ww-model-selection

This repository includes code for simulating the effect of damage to wastewater network components on the functionality of the network. Specifically, two methods are used for assessing the functionality, with two corresponding scripts here:
* flow.py: uses flow analysis (specifically SWMM) for functionality assessment
* connectivity.py: uses connectivity analysis for functionality assessment
In both cases, simulations continue until [convergence criteria] is met. 

# Specifying the Input Model

# Specifying the Damage Probabilities

# Outputs
Both scripts produce two outputs that define the network:
* ./outputs/nodes/nodes.shp: nodes shapefile that includes the following attributes
* ./outputs/edges/edges.shp: edges shapefile that includes the following attributes

# Reference
 Dunton, A. (2023). Model selection for risk analysis of wastewater networks. https://doi.org/10.48550/arXiv.2312.06623
