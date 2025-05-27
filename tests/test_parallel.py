from ferm import multiprocessing 
import scipy.sparse as sp

############LOAD DATA#########################################################
num_particles = 100
sigma = 1.0
path_niche_array = "../data/test_data/array_of_niche_to_pop_outest.npy"
path_x = "../data/test_data/x_pop.npy"
path_y = "../data/test_data/y_pop.npy"
path_pop = "../data/test_data/pop_test.tif"



# =============================================================================
P_final = multiprocessing(num_particles, sigma, path_niche_array,path_x,path_y,path_pop)
P_final = P_final.tocsr()
P_final = P_final / num_particles
sp.save_npz("test_sparse_mobiltiy_mat",P_final)
# =============================================================================
# for i in range(len(P_final)):
#     print(P_final[i].nnz)
# =============================================================================


