
clc;
clear all;

for k = 2003:1:2003
    clear tifImage;
%     inputFilename = sprintf('climate/sea_surface_temperature/sst_%d_raw.tif',k);
%     outputFilename = sprintf('climate/sea_surface_temperature/sst_%d.jpeg',k);
    inputFilename = sprintf('climate/ocean_acidification/slr_%d_impact.tif',k);
    outputFilename = sprintf('climate/ocean_acidification/slr_10_%d.jpeg',k);
    tifImage = imread(inputFilename);
    % 
    ii =0;
    jj =0;
        for i = 1:10:size(tifImage,1)
            ii = ii+1;
            for j = 1:10:size(tifImage,2)
                jj = jj+1;
                if (tifImage(i,j) < 0)
                     mage(ii,jj) = 0.00001;
                  else
                     mage(ii,jj) = tifImage(i,j);
                end    
            end
        end
%     rgb(:,:,1) = double(tifImage);
%     rgb(:,:,2) = zeros(size(tifImage,1),size(tifImage,2));
%     rgb(:,:,3) = zeros(size(tifImage,1),size(tifImage,2));
%         %
     imwrite(mage, outputFilename);
end