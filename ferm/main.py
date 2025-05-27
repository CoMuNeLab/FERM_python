import numpy as np
from geopy.distance import geodesic
from geokernels.distance import geodist
import scipy.sparse as sp
from tqdm import tqdm
import rioxarray
from multiprocessing import Pool


from ferm.utils import sample_max_distribution, gaussian_distribution_max

#### Compute niche which gives the mean of the gaussian algorithm

def distance_matrix(mask):
    mask_x = mask[0]
    mask_y = mask[1]
    distance_matrix = np.zeros([len(mask_x),len(mask_y)])
    for i in tqdm(range(len(mask_x))):
        for j in range(i,len(mask_x)):
            x,y = mask_x[i],mask_y[i]
            a,b = mask_x[j],mask_y[j]
            distance_matrix[i][j] = geodesic((x,y),(a,b))
    return distance_matrix



def parse_lat_lon(mask,x_pop,y_pop):
    mask_x = mask[0]
    mask_y = mask[1]
    coord_x = x_pop[mask_x]
    coord_y = y_pop[mask_y]
    points = np.array([(co_y,co_x) for co_x,co_y in zip(coord_x,coord_y)])
    return mask_x,mask_y,points


def wrap_geodist(a,b):
    return geodist(a,b,metric='km')


def FERM(
    nb_particules: int,
    sigma: float,
    path_niche_array: str,
    path_x: str,
    path_y: str,
    path_pop: str
        ) -> None:
    
    df_pop = rioxarray.open_rasterio(path_pop)
    # array of niche index corresponds to index in x_pop,y_pop
    array_niche = np.load(path_niche_array)
    x_pop = np.load(path_x)  # longitude coordiantes
    y_pop = np.load(path_y) # latitude coordiantes
    #print(len(x_pop),len(y_pop),array_niche.shape)
    array_niche = precise_the_mask(-40, 47, 4, x_pop, y_pop, array_niche)
    mask = np.where(array_niche !=0)
    mask_x,mask_y,points = parse_lat_lon(mask,x_pop,y_pop)
    P_final = sp.lil_matrix(((len(mask_x),len(mask_x)))) # the array of mobility from current point to every_destination
    for i in tqdm(range(len(mask_x))): #add here on all networks
        #print(i)
        p1 = np.array(points[i])
        x_current = mask[0][i]
        y_current = mask[1][i]
        pop_i = int(df_pop.sel(x=p1[1],y=p1[0]).values[0])
        if pop_i <1: #test if there is  population because if not move on
            continue
            #print("No population")
        else:
            distance_from_xy =  np.array([geodist(p1,points[j],metric='km') for j in range(len(points))])
            index_sort = np.argsort(distance_from_xy) #undoable for 1 million points
            index_sort = index_sort[1:] #the first element is the node itself so thrash
            for pt in range(nb_particules): # each particule coming from i
                #print(pt)
                mu = array_niche[x_current][y_current]
                absorption_i = gaussian_distribution_max(sigma,mu,pop_i)
                #print(absorption_i)
                for index in index_sort: #index of the point under scurtiny
                    p_destination = np.array(points[index])
                    pop_j = int(df_pop.sel(x=p_destination[1],y=p_destination[0]).values[0])
                    if pop_j >=1:
                        #print(len(index_sort))
                        x_destination = mask[0][index]
                        y_destination = mask[1][index]
                        mu_j =  array_niche[x_destination][y_destination]
                        #print(mu_j," lol ",pop_j)
                        absorbance_j = gaussian_distribution_max(sigma,mu_j,pop_j)
                        if absorbance_j > absorption_i:
                            P_final[i,index] += 1 #a lot of zeros in there, might ty to use lil matrix 
                            break # if absorption then break the for loop
                    else: #suppress index where there is nobody to gain time for the next particule
                         index_to_remove = np.where(index_sort == index)
                         index_sort = np.delete(index_sort,index_to_remove[0][0])
    P_final = P_final.tocsr()
    P_final = P_final / nb_particules
    sp.save_npz("pop=test_sparse_mobiltiy_mat",P_final)
    
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
            mu = array_niche[x_current][y_current]
            absorption_i = gaussian_distribution_max(sigma,mu,pop_i)
            #print(absorption_i)
            for index in index_sort: #index of the point under scurtiny
                p_destination = np.array(points[index])
                pop_j = int(df_pop.sel(x=p_destination[1],y=p_destination[0]).values[0])
                if pop_j >=1:
                    #print(len(index_sort))
                    x_destination = mask[0][index]
                    y_destination = mask[1][index]
                    mu_j =  array_niche[x_destination][y_destination]
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
    #print(args)
    with Pool(12,initializer,()) as pool:
        result = pool.starmap(FERM_multiprocessing,args)
    for i in range(len(result)):
        #print(result[i].nnz)
        P_final[i] = result[i]
    return P_final

def precise_the_mask(xmin, xmax, ymin, x_pop, y_pop, array_niche):
    x_cut = np.where(x_pop <=xmax,x_pop,abs(x_pop*0))
    x_cut = np.where(x_cut >=xmin,x_cut,abs(x_cut*0))
    y_cut = np.where(y_pop >=ymin,y_pop,abs(y_pop*0))
    print(np.where(y_cut==0))
    for x in np.where(x_cut==0)[0]:
        array_niche[x] = np.zeros(len(y_pop))
    
    for y in np.where(y_cut==0)[0]:
        array_niche[:,y] = np.zeros(len(x_pop))
    return array_niche

############LOAD DATA#########################################################
# nb_particules = 100
# sigma = 1
# path_niche_array = "../data/test_data/array_of_niche_to_pop_outest.npy"
# path_x = "../data/test_data/x_pop.npy"
# path_y = "../data/test_data/y_pop.npy"
# path_pop = "../data/test_data/pop_test.tif"

# import warnings
# warnings.filterwarnings("ignore")


    
    
#######################RUN#############################
# FERM(path_niche_array,path_x,path_y,path_pop)
# =============================================================================
# P_final = multiprocessing(path_niche_array,path_x,path_y,path_pop)
# P_final = P_final.tocsr()
# P_final = P_final / nb_particules
# sp.save_npz("test_sparse_mobiltiy_mat",P_final)
# =============================================================================
# =============================================================================
# for i in range(len(P_final)):
#     print(P_final[i].nnz)
# =============================================================================
# =============================================================================
# df_pop = rioxarray.open_rasterio(path_pop)
# array_niche = np.load(path_niche_array) #array of niche index corresponds to index in x_pop,y_pop
# x_pop = np.load(path_x) # longitude coordiantes
# y_pop = np.load(path_y) # latitude coordiantes
# mask = np.where(array_niche>0)
# mask_x,mask_y,points = parse_lat_lon(mask,x_pop,y_pop)
# P_final = sp.lil_matrix(((len(mask_x),len(mask_x)))) # the array of mobility from current point to every_destination
# for i in range(10000): #add here on all networks
#     print(i)
#     p1 = np.array(points[i])
#     x_current = mask[0][i]
#     y_current = mask[1][i]
#     pop_i = int(df_pop.sel(x=p1[1],y=p1[0]).values[0])
#     if pop_i <1: #test if there is  population because if not move on
#         continue
#         #print("No population")
#     else:
#         distance_from_xy =  np.array([geodist(p1,points[j],metric='km') for j in range(len(points))])
#         index_sort = np.argsort(distance_from_xy) #undoable for 1 million points
#         index_sort = index_sort[1:] #the first element is the node itself so thrash
#         for _ in range(nb_particules): # each particule coming from i
#             mu = array_niche[x_current][y_current]
#             absorption_i = gaussian_distribution_max(sigma,mu,pop_i)
#             #print(absorption_i)
#             for index in index_sort: #index of the point under scurtiny
#                 p_destination = np.array(points[index])
#                 pop_j = int(df_pop.sel(x=p_destination[1],y=p_destination[0]).values[0])
#                 if pop_j >=1:
#                     #print(len(index_sort))
#                     x_destination = mask[0][index]
#                     y_destination = mask[1][index]
#                     mu_j =  array_niche[x_destination][y_destination]
#                     #print(mu_j," lol ",pop_j)
#                     absorbance_j = gaussian_distribution_max(sigma,mu_j,pop_j)
#                     if absorbance_j > absorption_i:
#                         P_final[i,index] += 1 #a lot of zeros in there, might ty to use lil matrix 
#                         break # if absorption then break the for loop
#                 else: #suppress index where there is nobody to gain time for the next particule
#                      index_to_remove = np.where(index_sort == index)
#                      index_sort = np.delete(index_sort,index_to_remove[0][0])
# P_final = P_final.tocsr()
# P_final = P_final / nb_particules
# sp.save_npz("test_sparse_mobiltiy_mat",P_final)
# =============================================================================

   
