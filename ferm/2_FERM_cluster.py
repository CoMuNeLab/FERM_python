# /// script
# dependencies = [
#  "numpy == 2.1.3",
#  "rasterio == 1.4.3",
#  "rioxarray == 0.18.2",
#  "geokernels == 0.2.2",
#  "scipy == 1.15.2",
#  "ARSpy == 0.4",
# ]
# ///
import numpy as np
from geokernels.distance import geodist
import scipy.sparse as sp
import rioxarray
from multiprocessing import Pool
import os

from arspy.ars import adaptive_rejection_sampling

from utils import sample_max_distribution, gaussian_distribution_max

def test(sigma,n1,n2):
    s1,s2 =0,0
    for i in range(1000):
        z1= gaussian_distribution_max(sigma,0,n1)
        z2= gaussian_distribution_max(sigma,0.5,n2)
        if z1 > z2:
            s1 +=1
        else:
            s2+=1
    return s1,s2
            
        

def expectation(mu,sigma,n):
    lower_bound = mu+ sigma*np.sqrt(np.log(n))*(1/(np.pi*np.log(2)))
    upper_bound = mu+ sigma*np.sqrt(2*np.log(n))
    return lower_bound,upper_bound
    

#### Compute niche which gives the mean of the gaussian algorithm

def parse_lat_lon(mask,x_pop,y_pop):
    mask_x = mask[0]
    mask_y = mask[1]
    coord_x = x_pop[mask_x]
    coord_y = y_pop[mask_y]
    points = np.array([(co_y,co_x) for co_x,co_y in zip(coord_x,coord_y)])
    return mask_x,mask_y,points


def wrap_geodist(a,b):
    return geodist(a,b,metric='km')

    
def FERM_multiprocessing(path_pop,i,nb_particules,sigma):
    df_pop = rioxarray.open_rasterio(path_pop)
    p1 = np.array(points[i])
    x_current = mask[0][i]
    y_current = mask[1][i]
    pop_i = int(df_pop.sel(x=p1[1],y=p1[0]).values[0])
    P_process = sp.lil_matrix(((1,len(mask_x))))
    if pop_i <1: #test if there is  population because if not move on
        return P_process
    else:
        distance_from_xy =  np.array([geodist(p1,points[j],metric='km') for j in range(len(points))])
        index_sort = np.argsort(distance_from_xy) #undoable for 1 million points
        index_sort = index_sort[1:] #the first element is the node itself so thrash
        for _ in range(nb_particules): # each particule coming from i
            mu = array_niche[x_current][y_current]+1e-6
            absorption_i = gaussian_distribution_max(sigma,mu,pop_i)
            #print(absorption_i)
            for index in index_sort: #index of the point under scurtiny
                p_destination = np.array(points[index])
                pop_j = int(df_pop.sel(x=p_destination[1],y=p_destination[0]).values[0])
                if pop_j >=1:
                    #print(len(index_sort))
                    x_destination = mask[0][index]
                    y_destination = mask[1][index]
                    mu_j =  array_niche[x_destination][y_destination]+1e-6
                    #print(mu_j," lol ",pop_j)
                    absorbance_j = gaussian_distribution_max(sigma,mu_j,pop_j)
                    if absorbance_j > absorption_i:
                        P_process[0,index] += 1 #a lot of zeros in there, might ty to use lil matrix 
                        break # if absorption then break the for loop
                else: #suppress index where there is nobody to gain time for the next particule
                     index_to_remove = np.where(index_sort == index)
                     index_sort = np.delete(index_sort,index_to_remove[0][0])
    return P_process


def initializer():
    global array_niche
    global x_pop
    global y_pop
    global mask
    global mask_x
    global mask_y
    global points
    #global df_pop
    
def multiprocessing(path_niche_array,path_x,path_y,path_pop):
    #set_start_method('fork')
    global array_niche
    global x_pop
    global y_pop
    global mask
    global mask_x
    global mask_y
    global points
    #global df_pop
    #df_pop = rioxarray.open_rasterio(path_pop)
    array_niche = np.load(path_niche_array) #array of niche index corresponds to index in x_pop,y_pop
    x_pop = np.load(path_x) # longitude coordiantes
    y_pop = np.load(path_y) # latitude coordiantes
    mask = np.where(array_niche != 0) #eliminate place with no niche
    mask_x,mask_y,points = parse_lat_lon(mask,x_pop,y_pop)
    P_final = sp.lil_matrix(((len(mask_x),len(mask_x))))
    args = [[path_pop,i,nb_particules,sigma] for i in range(len(mask_x))]
    #FERM_multiprocessing(100)
    #print(args
    with Pool(os.cpu_count(),initializer,()) as pool:
        print("Number of CPU", os.cpu_count())
        result = pool.starmap(FERM_multiprocessing,args,chunksize=1)
    for i in range(len(result)):
        #print(result[i].nnz)
        P_final[i] = result[i]
    return P_final
    

############LOAD DATA#########################################################
nb_particules = 500
sigma = 1
main_path = "2070/"
path_niche_array = main_path+"array_of_niche_to_pop.npy"
path_x = main_path+"x_pop.npy"
path_y = main_path+"y_pop.npy"
path_pop = main_path+"pop=cluster_test.tif"
save_mobility = "mobility_sigma=1_cunksize=1"

import warnings
warnings.filterwarnings("ignore")
#######################RUN#############################

P_final = multiprocessing(path_niche_array,path_x,path_y,path_pop)
P_final = P_final.tocsr()
P_final = P_final / nb_particules
sp.save_npz(save_mobility,P_final)


   
