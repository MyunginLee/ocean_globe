
clc;
clear all;

for k = 2004:1:2004
    clear tifImage;
%     inputFilename = sprintf('climate/sea_surface_temperature/sst_%d_raw.tif',k);
%     outputFilename = sprintf('climate/sea_surface_temperature/sst_%d.jpeg',k);
%     inputFilename = sprintf('climate/ocean_acidification/slr_%d_impact.tif',k);
%     outputFilename = sprintf('climate/ocean_acidification/slr_10_%d.jpeg',k);
    inputFilename = sprintf('../texture/sst/equi_sst_%d.jpeg',k);
    outputFilename = sprintf('../texture/sst/hsv_equi_sst_%d.jpeg',k);

    mage = imread(inputFilename);

    for i = 1:1:351
        for j = 1:1:864
            c(i,j,1) = mage(i,j,1);
            c(i,j,2) = 0;
            c(i,j,3) = 0;
            d(i,j) = mage(i,j,1);
        end
    end


    imwrite(c, outputFilename);


end