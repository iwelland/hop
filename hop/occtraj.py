import hop.trajectory
import hop.graph
import numpy as np
import matplotlib.pyplot as plt

def occupancy_trajectory(hopgraph,trajectory):
    
    '''
    Returns a trajectory of the water network occupancies from a hop trajectory. Inputs are an unfiltered HoppingGraph and a  HoppingTrajectory object.
    '''

    n_nodes = len(hopgraph.graph.nodes())
    totaltime = trajectory.totaltime
    occ_traj = np.zeros((totaltime,n_nodes))
    
    for ts in trajectory:
        positions = ts.positions.T[0]
        positions_nonzero = positions[positions != 0].astype(int)
        occ_traj[ts.frame][positions_nonzero-1]=1
        #positions_zero = positions[positions == 0].astype(int)
        occ_traj[ts.frame][occ_traj[ts.frame]==0] = -1
        
    return occ_traj

def pair_correlations(occupancy_trajectory):

    occtraj=occupancy_trajectory
    totaltime = len(occtraj.T[0])
    norm = np.divide(1,totaltime,dtype=np.float32) 
    corr_matrix = norm*np.dot(occtraj.T,occtraj)
    return corr_matrix

def plot_timeseries(occupancy_trajectory,stride,separation,sites):

    index = 0
    for i in sites:
        traj = occupancy_trajectory.T[i]
        plt.plot(traj[:-1:stride]+separation*index,'p')
        index +=1

def path_clustering(occupancy_trajectory):
    
    dims = np.shape(occupancy_trajectory.T)
    means = np.empty((dims[0]))

    for i in xrange(dims[0]):
        occvect = occupancy_trajectory.T[i]
        means[i] = occvect.mean()

    def dist(mean1,mean2):
        return (mean1 - mean2)

    #dist_vec = np.vectorize(dist)
    dist_matrix = np.empty((dims[0],dims[0]))
    for i in xrange(dims[0]):
        for j in xrange(dims[0]):
            dist_matrix[i][j] = dist(means[i],means[j])

    return dist_matrix







    
