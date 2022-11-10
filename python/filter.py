import numpy as np
from math import sin, cos, pi, sqrt, asin
from PIL import Image
from scipy import interpolate
import cv2

Image.MAX_IMAGE_PIXELS = None

# years
for year in range(2003, 2014):
    print(year)
    path = '../texture/sst/test_median/test_0.5_sst_{}.png'.format(year) 
    img = Image.open(path)
    width, height = img.size
    img_data = np.array(img)

    blur = cv2.medianBlur(img_data, 5)
    im = Image.fromarray(blur)
    # rgb_img = Image.fromarray((np.dstack((img_data))).astype(np.uint8)) # recombing each color layer into an image
    path = '../texture/sst/test_median/sst_med_{}.png'.format(year) # https://map-projections.net/img/jpg/mollweide.jpg
    # img_data.save(path)
    im.save(path)
    # cv2.imwrite(path, proj)
    # rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/nutrient_pollution/image/equirectangular_nutrient_pollution_2003_impact.png')
    # rgb_img.save('/home/ben/Desktop/Projects/Sensorium/data/cumulative_human_impacts/climate/ocean_acidification/equirectangular_slr_2003.png')
    del img_data, blur, im
else:
    print("done")
