
clc;
clear all;

for k = 2003:1:2003
    clear tifImage;
%     inputFilename = sprintf('climate/sea_surface_temperature/sst_%d_raw.tif',k);
%     outputFilename = sprintf('climate/sea_surface_temperature/sst_%d.jpeg',k);
    inputFilename = sprintf('climate/ocean_acidification/slr_%d_impact.tif',k);
    outputFilename = sprintf('climate/ocean_acidification/slr_%d.csv',k);
    tifImage = imread(inputFilename);
    % 
        for i = 1:size(tifImage,1)
            for j = 1:size(tifImage,2)
                if (tifImage(i,j) < 0)
                     tifImage(i,j) = 0.00001;
%                  else
%                       tifImage(i,j) = tifImage(i,j);
                end    
            end
        end
        writematrix(tifImage, outputFilename)
        %
%     imwrite(tifImage, outputFilename);
end