"""
Author: Aaron Dunton (ajdunton@gmail.com)
"""

from swmm_api import SwmmInput, read_rpt_file
from swmm_api.input_file.macros.edit import delete_node, delete_link
from pyswmm import Simulation
import networkx as nx
import geopandas as gpd
import time
import random
import math
import sys, os
from contextlib import contextmanager
from shapely.geometry import LineString

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
        
    graph = nx.DiGraph()
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

def reduce_diams(model, pipe_damage):
    """
    Modifies some of the pipe diameters in model (SwmmInput, swmm_api), 
    according to the reduction factors in pipe_damage.

    Parameters
    ----------
    model : SwmmInput
        Input file as read by swmm-api.
    pipe_damage : dictionary
        Keys are conduits, and values are the corresponding reduction factors 
        for the diameters.

    Returns
    -------
    model : SwmmInput
        Modified input file, using swmm-api package.
    """
    
    for pipe, reduction in pipe_damage.items():
        model.XSECTIONS[str(pipe)].height = \
            reduction*model.XSECTIONS[str(pipe)].height

    return model

def sim(name):
    """
    Function that runs one simulation and saves results (performance measure) 
    for each node as column pm_+name in the nodes geodataframe. 

    Parameters
    ----------
    name : string
        Name of this simulation, pm_+name is column that is added to nodes 
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
    damaged_pipes = [] #[((from_node,to_node), fraction of flow cap), ...]
    
    for i in inp_undamaged.PUMPS.values():
        if random.random()<PF_PUMP:
            failed_pumps.append((int(i.from_node),int(i.to_node)))
        else:
            working_pumps.add((int(i.from_node),int(i.to_node)))
    for i in inp_undamaged.STORAGE.keys():
        if random.random()<PF_TANK:
            failed_tanks.append(int(i))
    for i in inp_undamaged.CONDUITS.values():
        ran_val = random.random()
        capacity_fractions = [i for i,_ in PF_PIPE]
        thresholds = [j for _,j in PF_PIPE]
        if ran_val > sum(thresholds[1:3]):
            if capacity_fractions[0]!=1 and capacity_fractions[0]!=0:
                damaged_pipes.append(((int(i.from_node),int(i.to_node)),
                                      capacity_fractions[0]))
            elif capacity_fractions[0]==0:
                failed_pipes.append((int(i.from_node),int(i.to_node)))
        elif ran_val > sum(thresholds[2:3]):
            if capacity_fractions[1]!=1 and capacity_fractions[1]!=0:
                damaged_pipes.append(((int(i.from_node),int(i.to_node)),
                                      capacity_fractions[1]))
            elif capacity_fractions[1]==0:
                failed_pipes.append((int(i.from_node),int(i.to_node)))
        elif ran_val > thresholds[3]:
            if capacity_fractions[2]!=1 and capacity_fractions[2]!=0:
                damaged_pipes.append(((int(i.from_node),int(i.to_node)),
                                      capacity_fractions[2]))
            elif capacity_fractions[2]==0:
                failed_pipes.append((int(i.from_node),int(i.to_node)))
        else: 
            if capacity_fractions[3]!=1 and capacity_fractions[3]!=0:
                damaged_pipes.append(((int(i.from_node),int(i.to_node)),
                                      capacity_fractions[3]))
            elif capacity_fractions[3]==0:
                failed_pipes.append((int(i.from_node),int(i.to_node)))
                
    # Modify, save, and run the input file    
    inp_damaged = inp_undamaged.copy()
    
    # Remove damaged pumps, storages, and all upstream items
    removed_nodes = set()
    removed_edges = set()
    edges_to_remove = [i for i in failed_pumps if i not in working_pumps] \
        + failed_pipes
    nodes_to_remove  = failed_tanks.copy()

    # Remove damaged pumps
    while len(edges_to_remove):
        cur_edge = edges_to_remove.pop(0)
        if cur_edge not in removed_edges:
            delete_link(inp_damaged,str(cur_edge))
            removed_edges.add(cur_edge)
            # Add upstream node to list of nodes to be removed
            nodes_to_remove.append(cur_edge[0])

    # Remove damaged storages
    while len(nodes_to_remove):  
        cur_node = nodes_to_remove.pop(0)
        if cur_node not in removed_nodes:
            delete_node(inp_damaged,str(cur_node))
            try: inp_damaged.DWF.pop((str(cur_node),'FLOW'))
            except KeyError: pass
            removed_nodes.add(cur_node)
            for k in g_undamaged.predecessors(cur_node):
                # Add upstream nodes to the list of nodes to be removed
                nodes_to_remove.append(k)
    
    if len(inp_damaged.CONDUITS)==0:
        nodes_gdf['pm_'+name] = 0
        nodes_gdf.loc[nodes_gdf.type=='outfall','pm'+name]=1
        return            
    
    # Reduce diamater of damaged pipes
    for pipe, reduction in damaged_pipes:
        if pipe[0] not in removed_nodes and pipe[1] not in removed_nodes:
            inp_damaged.XSECTIONS[str(ij_to_e[pipe])].height = \
                reduction*inp_damaged.XSECTIONS[str(ij_to_e[pipe])].height

    # Write modified SWMM file
    with suppress_stdout():
        inp_damaged.write_file('damaged.inp')
        sim_dam = Simulation('damaged.inp')
        sim_dam.execute()
    
    # Evaluate the performance measure from the SWMM output and save as column 
    # in nodes_gdf as column = name
        rpt_damaged = read_rpt_file('damaged.rpt')
    
    # Create damaged graph and add the max/full depth attribute to each of the
    # conduit edges
    g_damaged = g_undamaged.copy()
    g_damaged.remove_edges_from(removed_edges)
    g_damaged.remove_nodes_from(removed_nodes)
    
    for e, perc in zip(rpt_damaged.link_flow_summary.reset_index().Link, 
                       rpt_damaged.link_flow_summary['Max/_Full_Depth']):
        
        if not math.isnan(perc):
            g_damaged[e_to_ij[int(e)][0]][e_to_ij[int(e)][1]]['pflow'] = perc
        else: 
            g_damaged[e_to_ij[int(e)][0]][e_to_ij[int(e)][1]]['pflow'] = None
    
    # Performance measure for the nodes that have been removed (i.e., 
    # disconnected nodes) is 0
    for i in removed_nodes:
        nodes_gdf.loc[nodes_gdf.node==i,'pm_'+name] = 0
    
    # For other nodes, evaluate the performance measure based on the percentage
    # flow values along the path in the graph
    for i in g_damaged.nodes():
        pipe_count = 0
        crit_pipe_count = 0
        cur_node = i
        while len([j for j in g_damaged.successors(cur_node)])!=0:
            next_node = next(g_damaged.successors(cur_node))
            if g_damaged[cur_node][next_node]['pflow'] is not None:
                pipe_count = pipe_count + 1
                if g_damaged[cur_node][next_node]['pflow']>0.8:
                    crit_pipe_count = crit_pipe_count + 1
            cur_node = next_node
        if pipe_count != 0:      
            nodes_gdf.loc[nodes_gdf.node==i,'pm_'+name] = \
                1-crit_pipe_count/pipe_count
            
    os.remove('./damaged.inp')
    os.remove('./damaged.rpt')
    os.remove('./damaged.out')
    
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
        if i%10==0: print(f'Starting simulation {i}. Elapsed time {time.time()-start_time}.')
        sim(str(i))
        cols.append('pm_'+str(i))
    return cols

def take_avg(row,cols):
    """
    Evaluates the new average for the performance measure given a list of new 
    simulations to incorporate. 

    Parameters
    ----------
    row : row of a geodataframe
        Row must have N, pm_avg, and columns for each of the new simulations 
        listed in cols.
    cols : list of strings
        List of new simulation results that we are incorporatign into the 
        average.

    Returns
    -------
    float
        New average performance measure for this row.

    """
    total_val=row.N*row.pm_avg
    for i in cols:
        total_val = total_val+row[i]
    return total_val/(len(cols)+row.N)

def cov(row):
    if row.pm_avg==0:
        val = 0
    else:
        val = math.sqrt((1-row.pm_avg)/(row.N*row.pm_avg))
    return val

def avg_cov():
    cov_series = nodes_gdf.apply(cov,axis=1)
    return cov_series.mean()

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout


# %% CONSTANTS
# Define damage probabilities for the pumps, pipes, and tanks
PF_PUMP = 0.1
PF_TANK = 0.1

# Fraction of flow capacity for DS : probability of DS
# Currently, the sim function assumes 4 damage states
PF_PIPE = [(1,0.85), 
           (0.5,0.1), 
           (0.1,0.04), 
           (0,0.01)] 

# Input file location
INP_PATH = './inputs/arterial_dynwave/arterial_dynwave.inp'

# %% UNDAMAGED MODEL
# Open the input file and create a graph
with suppress_stdout():
    inp_undamaged = SwmmInput.read_file(INP_PATH)
outfall_ind, g_undamaged = graph_from_SwmmInput(inp_undamaged)
e_to_ij = dict()
for edge in inp_undamaged.CONDUITS.values():
    e_to_ij[int(edge.name)] = (int(edge.from_node),int(edge.to_node))
for edge in inp_undamaged.PUMPS.values():
    e_to_ij[int(edge.name)] = (int(edge.from_node),int(edge.to_node))
ij_to_e = {val:key for key,val in e_to_ij.items()}

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
nodes_gdf['pm_avg'] = 0
cur_N = 0
new_cols = []

while avg_cov()>0.1 or cur_N<100:
    if cur_N > 0: print(f'Current average cov is {avg_cov()}.')
    if cur_N > 0: print(f'Current elapsed time is {time.time()-start_time}.')
    nodes_gdf = nodes_gdf.drop(new_cols, axis=1)
    new_cols = run_100(cur_N)
    nodes_gdf['pm_avg'] = nodes_gdf.apply(take_avg, axis=1, cols=new_cols)
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