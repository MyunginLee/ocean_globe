import numpy as np
from math import sin, cos, pi, sqrt, asin
from PIL import Image
from scipy import interpolate
import cv2
from timeit import default_timer as timer
import multiprocessing as mp
from multiprocessing import Process, Value, Array
import os

Image.MAX_IMAGE_PIXELS = None

sqt2 = sqrt(2)
def inv_mollweide(x, y):
    """Mollweide inverse tranformation"""
    
    theta = asin(y/(sqt2)) # auxiliary angle
    
    lat = asin((2*theta+sin(2*theta))/pi) # latitude
    lon = (pi*x)/(2*sqt2*cos(theta)) # longitude
    
    return lat, lon

def moolweide2geodetic(data, q):
    """Applies Mollweide inverse tranformation on each cell of the grid"""
    
    n_lat, n_lon = data.shape # retrieving number of rows and columns
    
    x = np.linspace(-2*sqt2, 2*sqt2, n_lon) # x definition domain
    y = np.linspace(-1*sqt2, sqt2, n_lat) # y definition domain
    
    lon = np.linspace(-pi, pi, n_lon) # longitude definition domain
    lat = np.linspace(-pi/2, pi/2, n_lat) # latitude definition domain
    
    data_proj = np.zeros((n_lat, n_lon)) # creating new array
    
    r = n_lon/2
    h = n_lat/2
    for i in range(n_lat):
        for j in range(n_lon):
            # j_lat = int(r + (j-r) * (r-(2/ pi * sqrt(2*h*h-(i-h)*(i-h)))) / r)
            p = ( r * 2 * sqt2 / pi * sqrt( 1 - (i-h)*(i-h) / (2*r*r)) )
            j_lat = int(r + (j-r) * p/ r)
            if data[i, j_lat] < 0:
                data[i, j_lat] = -0.01
            else:
                data[i, j_lat] = data[i, j_lat] *255
            data_proj[i, j] = data[i, j_lat] # EDIT



    data_proj = np.ma.masked_invalid(data_proj) # EDIT
    
    LON, LAT = np.meshgrid(lon, lat) # EDIT
    LON1 = LON[~data_proj.mask] # EDIT
    LAT1 = LAT[~data_proj.mask] # EDIT
    data_proj_new = data_proj[~data_proj.mask] # EDIT
    
    GD = interpolate.griddata((LON1, LAT1), data_proj_new.ravel(),
                              (LON, LAT), method='linear') # EDIT
    # GD[np.isnan(GD)] = 0 # EDIT
    
    # return GD # EDIT
    q.put(GD)

if __name__ == '__main__':
    print('pid of main:', os.getpid())
    q = mp.Queue()
    # years
    # for year in range(2003, 2008):
    for year in range(2003, 2004):
        print(year)
        # path = 'C:/Users/allos/Desktop/db/organic_chem/organic_chemical_pollution_{}_impact.tif'.format(year)
        path = '/Users/ben/Desktop/projects/sensorium/chi_impact/sst/sst_{}_impact.tif'.format(year)
        out_path = '/Users/ben/Desktop/projects/sensorium/chi_impact/sst/sst_5_{}_impact.png'.format(year)

        start = timer()
        # Read image
        img = Image.open(path)

        img_data = np.array(img)

        # trim_img = img_data[0:18120,:]
        width, height = img_data.shape

        img = cv2.resize(img_data, dsize=(int(height/20),int(height/40)), interpolation=cv2.INTER_CUBIC)
        # tmp = moolweide2geodetic(img) # applying Mollweide inverse tranformation on the red grid
        # proj = np.array(tmp, dtype='uint8')

        p1 = Process(target=moolweide2geodetic, args = (img, q,))
        p1.start()
        tmp = q.get()
        p1.join()
        proj = np.array(tmp, dtype='uint8')

        im = Image.fromarray(proj)
        # path = 'C:/Users/allos/Desktop/db/organic_chem/oc_10_{}_impact.png'.format(year)
        # img_data.save(path)
        im.save(out_path)
        del img_data, im, tmp, proj
        print("time:", timer()-start)
    else:
        print("done")
