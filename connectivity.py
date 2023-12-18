"""
Author: Aaron Dunton (ajdunton@gmail.com)
"""

import networkx as nx
from swmm_api import SwmmInput
import random
import geopandas as gpd
from shapely.geometry import LineString
import math
import time

# %% FUNCTIONS
def graph_from_SwmmInput(sw_in):
    """
    Creates a NetworkX DiGraph corresponding the SWMM network, including all 
    nodal components (junctions, outfalls, storage) and all edge components 
    (conduits, pumps, orifices).

    Parameters
    ----------
    sw_in : SwmmInput 
        SWMM file as opened by swmm_api.

    Returns
    -------
    outfall : integer
        Node number for the outfall.
    graph : NetworkX DiGraph
        Directed graph representing the input SWMM network.
    """
    node_set = set()
    for i in sw_in.JUNCTIONS.keys():
        node_set.add(int(i))
    for i in sw_in.OUTFALLS.keys():
        node_set.add(int(i))
        outfall = int(i)
    for i in sw_in.STORAGE.keys():
        node_set.add(int(i))
        
    edge_set = set()
    for i in sw_in.CONDUITS.values():
        edge_set.add((int(i.from_node),int(i.to_node)))
    for i in sw_in.PUMPS.values():
        edge_set.add((int(i.from_node),int(i.to_node))) 
    if sw_in.ORIFICES is not None:
        for i in sw_in.ORIFICES.values():
            edge_set.add((int(i.from_node),int(i.to_node)))
        
    graph = nx.Graph()
    graph.add_nodes_from(node_set)
    graph.add_edges_from(edge_set)
    
    return outfall, graph

def add_type_to_nodes(row):
    if str(row.node) in inp_undamaged.STORAGE.keys():
        type = 'tank'
    elif str(row.node) in inp_undamaged.JUNCTIONS.keys():
        type = 'junction'
    elif str(row.node) in inp_undamaged.OUTFALLS.keys():
        type = 'outfall'
    return type

def sim(name):
    """
    Function that runs one simulation and saves results (indicator of connected
    or disconnected) for each node as column con_+name in the nodes 
    geodataframe. 

    Parameters
    ----------
    name : string
        Name of this simulation, con_+name is column that is added to nodes 
        gdf.

    Returns
    -------
    None.

    """
    
    # Component Damage
    failed_pumps = []
    working_pumps = set()
    failed_tanks = []
    failed_pipes = []
    
    for i in inp_undamaged.PUMPS.values():
        if random.random()<PF_PUMP:
            failed_pumps.append((int(i.from_node),int(i.to_node)))
        else:
            working_pumps.add((int(i.from_node),int(i.to_node)))
    final_failed_pumps = [i for i in failed_pumps if i not in working_pumps]
    for i in inp_undamaged.STORAGE.keys():
        if random.random()<PF_TANK:
            failed_tanks.append(int(i))
    for i in inp_undamaged.CONDUITS.values():
        if random.random()<PF_PIPE:
            failed_pipes.append((int(i.from_node),int(i.to_node)))
            
    # Connectivity Analysis
    g_damaged = g_undamaged.copy()
    g_damaged.remove_edges_from(final_failed_pumps)
    g_damaged.remove_edges_from(failed_pipes)
    g_damaged.remove_nodes_from(failed_tanks)

    if not nx.is_connected(g_damaged):
        nodes_gdf['con_'+name] = 0
        for comp in nx.connected_components(g_damaged):
            if outfall_ind in comp:
                for i in comp:
                    nodes_gdf.loc[nodes_gdf.node==i,'con_'+name] = 1
    else: nodes_gdf['con_'+name] = 1
    return

def run_100(start):
    """
    Runs 100 simulations, calling sim for each simulation. Keeps track of the 
    simulation number to tell sim what to call the new column and to return a 
    list of the new columns. 

    Parameters
    ----------
    start : int
        Simulation nunber for the first simulation in this batch of 100.

    Returns
    -------
    cols : list
        List of strings, names of new columns added to nodes GDF (i.e., one 
        column for each simulation).

    """
    cols = []
    print(f'Running simulations {start}-{start+99}.')
    for i in range(start,start+100):
        sim(str(i))
        cols.append('con_'+str(i))
    return cols

def take_avg(row,cols):
    """
    Evaluates the new average for the performance measure given a list of new 
    simulations to incorporate. 

    Parameters
    ----------
    row : row of a geodataframe
        Row must have N, con_avg, and columns for each of the new simulations 
        listed in cols.
    cols : list of strings
        List of new simulation results that we are incorporatign into the 
        average.

    Returns
    -------
    float
        New average performance measure for this row.

    """
    total_val=row.N*row.con_avg
    for i in cols:
        total_val = total_val+row[i]
    return total_val/(len(cols)+row.N)

def cov(row):
    if row.con_avg==0:
        val = 0
    else:
        val = math.sqrt((1-row.con_avg)/(row.N*row.con_avg))
    return val

def avg_cov():
    cov_series = nodes_gdf.apply(cov,axis=1)
    return cov_series.mean()

# %% CONSTANTS
# Define damage probabilities for the pumps, pipes, and tanks
PF_PUMP = 0.1
PF_TANK = 0.1
PF_PIPE = 0.01 

# Input file location
INP_PATH = './inputs/arterial_dynwave/arterial_dynwave.inp'

# %% UNDAMAGED MODEL
# Open the input file and create a graph
inp_undamaged = SwmmInput.read_file(INP_PATH)
outfall_ind, g_undamaged = graph_from_SwmmInput(inp_undamaged)

# Geodataframe of nodes, include DWF columns
dwf_df = inp_undamaged.DWF.get_dataframe()
nodes_gdf = \
    gpd.GeoDataFrame(inp_undamaged.COORDINATES.get_geo_series(crs=2913))
nodes_gdf = nodes_gdf.merge(dwf_df, 
                            left_index=True, 
                            how='left',
                            right_on='node').reset_index().astype({'node':int})
nodes_gdf = nodes_gdf.drop('index',axis=1)
nodes_gdf['type'] = nodes_gdf.apply(add_type_to_nodes,axis=1)

# %% SIMULATIONS
start_time = time.time()

nodes_gdf['N'] = 0
nodes_gdf['con_avg'] = 0
cur_N = 0
new_cols = []

while avg_cov()>0.1 or cur_N<100:
    if cur_N > 0: print(f'Current average cov is {avg_cov()}.')
    nodes_gdf = nodes_gdf.drop(new_cols, axis=1)
    new_cols = run_100(cur_N)
    nodes_gdf['con_avg'] = nodes_gdf.apply(take_avg, axis=1, cols=new_cols)
    cur_N = cur_N + 100
    nodes_gdf['N'] = cur_N
    
end_time = time.time()

print(f'Convergence criteria satisfied (cov = {avg_cov()}) after {cur_N} simulations.')
print(f'Elapsed time: {end_time-start_time}.')
    
# %% SAVE OUTPUT
# Save node gdf
nodes_gdf.to_file('./outputs/nodes')

# Save edge gdf
from_node = []
to_node = []
e_type = []
geom = []

for i,j in g_undamaged.edges():
    pt0 = nodes_gdf.loc[nodes_gdf.node==i].geometry.values[0]
    pt1 = nodes_gdf.loc[nodes_gdf.node==j].geometry.values[0]
    geom.append(LineString([pt0,pt1])) 
    from_node.append(i)
    to_node.append(j)
    if i in [int(k.from_node) for k in inp_undamaged.CONDUITS.values()]:
        e_type.append('pipe')
    elif i in [int(k.from_node) for k in inp_undamaged.PUMPS.values()]:
        e_type.append('pump')
    else:
        e_type.append('other')

edges_gdf = gpd.GeoDataFrame(data = {'from_n':from_node,
                                     'to_n':to_node,
                                     'type':e_type,
                                     'geometry':geom}).set_crs(2913)

edges_gdf.to_file('./outputs/edges')