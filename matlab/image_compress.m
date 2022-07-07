
clc;
clear all;

for k = 2003:1:2003
    clear tifImage;
%     inputFilename = sprintf('climate/sea_surface_temperature/sst_%d_raw.tif',k);
%     outputFilename = sprintf('climate/sea_surface_temperature/sst_%d.jpeg',k);
%     inputFilename = sprintf('climate/ocean_acidification/slr_%d_impact.tif',k);
%     outputFilename = sprintf('climate/ocean_acidification/slr_10_%d.jpeg',k);
    inputFilename = sprintf('nutrient_pollution/image/image_nutrient_pollution_%d_impact.jpeg',k);
    outputFilename = sprintf('nutrient_pollution/image/image_nutrient_pollution_10_%d.jpeg',k);

    mage = imread(inputFilename);

    [U1,S1,V1] = svdsketch(double(mage),1e-2);
    Anew1 = uint8(U1*S1*V1');


    imwrite(Anew1, outputFilename);


end