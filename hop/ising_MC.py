import numpy as np
import hop.graph
import random
import hop.walker

def return_sites(filtered_hopgraph,z_center):
    '''
    Returns a list of sites above z_center and sites below z_center. Should be added to a new module
    that is separate from hop.walker or hop.ising_MC.
    '''
    h=filtered_hopgraph
    if up_flux==True:
        top_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
        h.site_properties.center[site][2]>=z_center]
        top_sites=[site for site in top_sites_unfiltered if h.graph.has_edge(1,site)]
        bottom_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
        h.site_properties.center[site][2]<=z_center]
        bottom_sites=[site for site in bottom_sites_unfiltered if h.graph.has_edge(site,1)]
    else:
        top_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
        h.site_properties.center[site][2]<=z_center]
        top_sites=[site for site in top_sites_unfiltered if h.graph.has_edge(1,site)]
        bottom_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
        h.site_properties.center[site][2]>=z_center]
        bottom_sites=[site for site in bottom_sites_unfiltered if h.graph.has_edge(site,1)]
    return top_sites,bottom_sites

def bulk_rate_calculator(hopgraph,sites):
        '''
        This function computes a clock time which determines when particles should be added or removed from the bulk sites.
        It sums over all bulk edges of sites and takes the inverse of the sum of all rates.
        NOTE: This function was copied from hop.walker.flux_calculator, and should be incorporated
        into a separate module (possibily utilities.py) rather than copied
        '''

        h=hopgraph
        bulk_rate=0
        for site in h.graph[1]:
            if site in sites:
                bulk_rate+=h.graph[1][site]['k']

        bulk_time=1/bulk_rate
        return bulk_rate

def init_J_matrix(filtered_hopgraph,totaltime):

    f = filtered_hopgraph
    #try: 
    #    f.filtered_graph
    #except AttributeError:
    #    raise 

    n_nodes = len(f.filtered_graph.nodes())
    J_matrix = np.zeros((n_nodes,n_nodes))
    for n_i in xrange(n_nodes):
        node_i = f.filtered_graph.nodes()[n_i]
        for n_j in xrange(n_nodes):
            node_j = f.filtered_graph.nodes()[n_j]
            if f.filtered_graph.has_edge(node_i,node_j):
                J_matrix[n_i][n_j] = f.filtered_graph[node_i][node_j]['N']

    return J_matrix*np.float(1)/totaltime

def Hamiltonian(filtered_graph,bulk_rates,sites,J_matrix,n_rows,state):

    f = filtered_graph
    #g = hop.walker.graph_update(f)
    n_rows = J_matrix.shape[0]
    H = 0
    #sites = f.filtered_graph.nodes()
    #bulk_rate = bulk_rate_calculator(f,sites)
    #B = bulk_rate/len(sites)
    B=bulk_rates.sum()
    for i in xrange(n_rows):
        S_i = state[i]
        #H+=-B*S_i
        for j in xrange(n_rows):
            S_j = state[j]
            H += -0.25*J_matrix[i][j]*(S_j+1)*(S_i+1) #+ 0.5*bulk_rates[i]*(1-S_j)
#            H += -0.5*(J_matrix[i][j]*(1+S_i)-J_matrix[j][i]*(1+S_j))-B #+ 0.5*bulk_rates[i]*(1-S_j)
    return H

def Metro_MC(filtered_hopgraph,initial_configuration,J_matrix,n_steps,totaltime,stride=1):

    f = filtered_hopgraph
    n_nodes = len(initial_configuration)
    traj = np.zeros((n_steps/stride,n_nodes))
    stride_index = 0
    n_strides = 0
    sites = f.filtered_graph.nodes()
    bulk_rate = bulk_rate_calculator(f,sites)
    B = bulk_rate/len(sites)
    n_rows = J_matrix.shape[0]
    bulk_rates=np.zeros(len(f.graph.nodes()))
    for node in f.graph.nodes():
        if f.graph.has_edge(node,1):
                bulk_rates[node-1]=f.graph[node][1]['N']
    bulk_rates=bulk_rates*1/totaltime
    for step in xrange(n_steps):
        #print("This is the step")
        #print(step)
        #print("This is the stride_index")
        #print(stride_index)
        #print("This is n_strides")
        #print(n_strides)
        if step == 0:
            for index in xrange(n_nodes):
                if initial_configuration[index] == 0:
                    initial_configuration[index]=-1
            traj[0]=initial_configuration
            last_state=traj[0]
        else:
            site_index = random.choice(range(n_nodes))
            last_state = traj[n_strides-1]
            site = last_state[site_index]
            site_new = -1*site
            state_new = np.zeros(len(last_state))
            state_new[...]=last_state
            state_new[site_index]=site_new
            delta_H = Hamiltonian(f,bulk_rates,sites,J_matrix,n_rows,state_new)-Hamiltonian(f,bulk_rates,sites,J_matrix,n_rows,last_state)
            if delta_H < 0:
                last_state=state_new
                if stride_index==0:
                    traj[n_strides]=state_new
            elif delta_H >= 0:
                alpha = np.random.uniform()
                if alpha < np.exp(delta_H):
                    last_state=state_new
                    if stride_index==0:
                        traj[n_strides]=state_new
                else:
                    if stride_index==0:
                        traj[n_strides]=last_state
        stride_index+=1
        if stride_index == stride and step!=n_steps-1:
#            print(100*float(step)/float(n_steps) + (n_steps/stride)**-1)
            stride_index = 0
            n_strides+=1
    return traj
        
#def kinetic_MC(initial_configuration,max_num_particles):

 #   particles=np.zeros(max_num_particles)
  #  init = initial_configuration
    
