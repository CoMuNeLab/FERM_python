from ferm import FERM 


############LOAD DATA#########################################################
num_particles = 100
sigma = 1.0
path_niche_array = "../data/test_data/array_of_niche_to_pop_outest.npy"
path_x = "../data/test_data/x_pop.npy"
path_y = "../data/test_data/y_pop.npy"
path_pop = "../data/test_data/pop_test.tif"


##### RUN ######
FERM(num_particles, sigma, path_niche_array,path_x,path_y,path_pop)