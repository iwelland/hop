import numpy as np
import hop.graph
import random

def init_J_matrix(filtered_hopgraph,totaltime):

    f = filtered_hopgraph
    #try: 
    #    f.filtered_graph
    #except AttributeError:
    #    raise 

    n_nodes = len(f.filtered_graph.nodes())
    J_matrix = np.empty((n_nodes,n_nodes))
    for n_i in xrange(n_nodes):
        node_i = f.filtered_graph.nodes()[n_i]
        for n_j in xrange(n_nodes):
            node_j = f.filtered_graph.nodes()[n_j]
            if f.filtered_graph.has_edge(node_i,node_j):
                J_matrix[n_i][n_j] = f.filtered_graph[node_i][node_j]['N']

    return J_matrix*np.float(1)/totaltime

def Hamiltonian(J_matrix,state):

    n_rows = J_matrix.shape[0][0]
    H = 0
    for i in xrange(n_rows):
        S_i = state[i]
        for j in xrange(n_rows):
            S_j = state[j]
            H += J[i][j]*S_i*(1-S_j)
    return -H

def Metro_MC(filtered_hopgraph,initial_configuration,J_matrix,n_steps):

    f = filtered_hopgraph
    n_nodes = len(f.filtered_graph.nodes())
    traj = np.empty((n_steps,n_nodes))
    for step in xrange(n_steps):
        if step == 0:
            traj[0]=initial_configuration
        else:
            site_index = random.choice(n_nodes)
            state = traj[step]
            site = state[site_index]
            site_new = -1*site
            state_new = state
            state_new[site_index]=site_new
            delta_H = H(J_matrix,state_new)-H(J_matrix,state)
            alpha = np.random.uniform()
            if alpha < np.exp(-delta_H)/(1 + np.exp(-delta_H)):
                traj[step]=state_new
            else:
                traj[step]=state
    return traj
        
        
        
