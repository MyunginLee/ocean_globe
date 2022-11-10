import numpy as np
from math import sin, cos, pi, sqrt, asin
from PIL import Image
from scipy import interpolate
import cv2

Image.MAX_IMAGE_PIXELS = None

def inv_mollweide(x, y, lon_0, R):
    """Mollweide inverse tranformation"""
    
    theta = asin(y/(R*sqrt(2))) # auxiliary angle
    
    lat = asin((2*theta+sin(2*theta))/pi) # latitude
    lon = lon_0 + (pi*x)/(2*R*sqrt(2)*cos(theta)) # longitude
    
    return lat, lon

def moolweide2geodetic(data, lon_0, R):
    """Applies Mollweide inverse tranformation on each cell of the grid"""
    
    n_lat, n_lon = data.shape # retrieving number of rows and columns
    
    x = np.linspace(-2*R*sqrt(2), 2*R*sqrt(2), n_lon) # x definition domain
    y = np.linspace(-R*sqrt(2), R*sqrt(2), n_lat) # y definition domain
    
    lon = np.linspace(-pi, pi, n_lon) # longitude definition domain
    lat = np.linspace(-pi/2, pi/2, n_lat) # latitude definition domain
    
    data_proj = np.zeros((n_lat, n_lon)) # creating new array
    
    for i in range(n_lat):
        for j in range(n_lon):
            
            coord = inv_mollweide(x[j], y[i], lon_0, R) # applying Mollweide inverse tranformation on cell (i, j)
            
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
    
    return GD # EDIT

# years
for year in range(2003, 2014):
    print(year)
    # path = '../texture/sst/sst_2003.jpeg' # https://map-projections.net/img/jpg/mollweide.jpg
    path = '/Users/ben/Desktop/projects/sensorium/chi_impact/slr/slr_{}_impact.tif'.format(year)
    # path = '/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/image_nutrient_pollution_2003_impact.jpeg'
    # path = '/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/slr_2003.jpeg'
    # Read image
    img = Image.open(path)
    width, height = img.size
    img_data = np.array(img)
    # newsize = (int(width/30), int(height/30)) #decreasing the quality to speed up the testing process
    # newsize = (500, 250) #decreasing the quality to speed up the testing process
    # img = cv2.resize(img_data, dsize=(500,250), interpolation=cv2.INTER_CUBIC)
    img = cv2.resize(img_data, dsize=(1500,750), interpolation=cv2.INTER_CUBIC)
    lon_0, R = 0, 1
    tmp = moolweide2geodetic(img, lon_0, R) # applying Mollweide inverse tranformation on the red grid
    proj = np.array(tmp, dtype='uint8')
    # im = Image.fromarray(proj)
    # blur = cv2.GaussianBlur(proj, (3,3), 1,1)
    # blur = cv2.medianBlur(proj, 3)
    blur = proj
    im = Image.fromarray(blur)
    # rgb_img = Image.fromarray((np.dstack((img_data))).astype(np.uint8)) # recombing each color layer into an image
    path = '../texture/slr/slr_{}_impact.jpg'.format(year) # https://map-projections.net/img/jpg/mollweide.jpg
    # img_data.save(path)
    im.save(path)
    # cv2.imwrite(path, proj)
    # rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/equirectangular_nutrient_pollution_2003_impact.png')
    # rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/equirectangular_slr_2003.png')
    del img_data, blur, im, tmp
else:
    print("done")
