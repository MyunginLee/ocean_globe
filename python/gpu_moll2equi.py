import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import torch
import numpy as np
# import cupy as np
from numba import njit, cuda, vectorize, float64
from math import sin, cos, pi, sqrt, asin
from PIL import Image
from scipy import interpolate
from timeit import default_timer as timer

import cv2
torch.__version__
torch.cuda.is_available()
print(torch.cuda.is_available())
device = "cuda" if torch.cuda.is_available() else "cpu" 
print("Using {} device".format(device))

Image.MAX_IMAGE_PIXELS = None

@njit( parallel=True)
def inv_mollweide(x, y):
    """Mollweide inverse tranformation"""
    theta = asin(y/(sqrt(2))) # auxiliary angle
    
    lat = asin((2*theta+sin(2*theta))/pi) # latitude
    lon = (pi*x)/(2*sqrt(2)*cos(theta)) # longitude
    
    return lat, lon

# @jit(target_backend ="cuda")
# @cuda.jit
@njit( parallel=True)
def moolweide2geodetic(data):
    print("data: ",type(data))
    """Applies Mollweide inverse tranformation on each cell of the grid"""
    n_lat, n_lon = data.shape # retrieving number of rows and columns
    
    x = np.linspace(-2*sqrt(2), 2*sqrt(2), n_lon) # x definition domain
    y = np.linspace(-sqrt(2), sqrt(2), n_lat) # y definition domain
    
    lon = np.linspace(-pi, pi, n_lon) # longitude definition domain
    lat = np.linspace(-pi/2, pi/2, n_lat) # latitude definition domain
    
    data_proj = np.zeros((n_lat, n_lon)) # creating new array
    
    for i in range(n_lat):
        for j in range(n_lon):
            
            coord = inv_mollweide(x[j], y[i]) # applying Mollweide inverse tranformation on cell (i, j)
            
            i_lat = list(lat).index(min(lat, key=lambda x:abs(x-coord[0]))) # retrieving new row index
            i_lon = list(lon).index(min(lon, key=lambda x:abs(x-coord[1]))) # retrieving new colmn index
            
            data_proj[i_lat, i_lon] = data[i, j] # EDIT
    
    data_proj = np.ma.masked_invalid(data_proj) # EDIT
    
    LON, LAT = np.meshgrid(lon, lat) # EDIT
    LON1 = LON[~data_proj.mask] # EDIT
    LAT1 = LAT[~data_proj.mask] # EDIT
    data_proj_new = data_proj[~data_proj.mask] # EDIT
    
    GD = interpolate.griddata((LON1, LAT1), data_proj_new.ravel(),
                              (LON, LAT), method='linear') # EDIT
    
    GD[np.isnan(GD)] = 0 # EDIT
    print("GD: ", type(GD))

    return GD # EDIT

# years
start = timer()
for year in range(2013, 2014):
    print(year)
    # path = '../texture/sst/sst_2003.jpeg' # https://map-projections.net/img/jpg/mollweide.jpg
    # path = '/Users/ben/Desktop/projects/sensorium/chi_impact/direct_human/direct_human_{}_raw.tif'.format(year)
    # path = '/Users/ben/Desktop/projects/sensorium/chi_impact/nutrient_pollution/nutrient_pollution_{}_impact.tif'.format(year)
    path = 'C:/Users/hyo96/projects/sensorium/db/nutrient_pollution_{}_impact.tif'.format(year)
    # path = '/Users/ben/Desktop/projects/sensorium/chi_impact/sst/sst_{}_raw.tif'.format(year)
    # path = '/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/image_nutrient_pollution_2003_impact.jpeg'
    # path = '/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/slr_2003.jpeg'
    # Read image
    img = Image.open(path)
    width, height = img.size
    img_data = np.array(img)
    # newsize = (int(width/30), int(height/30)) #decreasing the quality to speed up the testing process
    #decreasing the quality to speed up the testing process 2:1
    # img = cv2.resize(img_data, dsize=(2574,1287), interpolation=cv2.INTER_CUBIC)
    # img = cv2.resize(img_data, dsize=(int(width/10),int(width/20)), interpolation=cv2.INTER_CUBIC)
    img = cv2.resize(img_data, dsize=(500,250), interpolation=cv2.INTER_CUBIC)
    tmp = moolweide2geodetic(img) # applying Mollweide inverse tranformation on the red grid
    proj = np.array(tmp, dtype='uint8')
    # im = Image.fromarray(proj)
    # blur = cv2.GaussianBlur(proj, (3,3), 1,1)
    # blur = cv2.medianBlur(proj, 5)
    blur = proj
    im = Image.fromarray(blur)
    # rgb_img = Image.fromarray((np.dstack((img_data))).astype(np.uint8)) # recombing each color layer into an image
    path = 'C:/Users/hyo96/projects/sensorium/db/np_{}.jpg'.format(year) # https://map-projections.net/img/jpg/mollweide.jpg
    # img_data.save(path)
    im.save(path)
    # cv2.imwrite(path, proj)
    # rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/equirectangular_nutrient_pollution_2003_impact.png')
    # rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/equirectangular_slr_2003.png')
    del img_data, blur, im, tmp
else:
    print("done")
print("with GPU:", timer()-start)
