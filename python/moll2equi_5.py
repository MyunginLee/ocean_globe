import numpy as np
from math import sin, cos, pi, sqrt, asin
from PIL import Image
from scipy import interpolate
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

path = '../texture/sst/sst_2012.jpeg' # https://map-projections.net/img/jpg/mollweide.jpg
# path = '/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/image_nutrient_pollution_2003_impact.jpeg'
# path = '/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/slr_2003.jpeg'

# Read image
img = Image.open(path)
width, height = img.size

# newsize = (int(width/30), int(height/30)) #decreasing the quality to speed up the testing process
newsize = (int(width/3), int(height/3)) #decreasing the quality to speed up the testing process
img = img.resize(newsize)

img_data = np.transpose(np.array(img), (2, 0, 1))
r = img_data[0, :, :]
g = img_data[1, :, :]
b = img_data[2, :, :]

lon_0, R = 0, 1

r_proj = moolweide2geodetic(r, lon_0, R) # applying Mollweide inverse tranformation on the red grid
g_proj = moolweide2geodetic(g, lon_0, R) # applying Mollweide inverse tranformation on the green grid
b_proj = moolweide2geodetic(b, lon_0, R) # applying Mollweide inverse tranformation on the blue grid


rgb_img = Image.fromarray((np.dstack((r_proj, g_proj, b_proj))).astype(np.uint8)) # recombing each color layer into an image
path = '../texture/sst/equi_sst_2012.jpeg' # https://map-projections.net/img/jpg/mollweide.jpg
rgb_img.save(path)

# rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/equirectangular_nutrient_pollution_2003_impact.png')
# rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/equirectangular_slr_2003.png')

del rgb_img, img_data, r_proj, g_proj, b_proj

