%% Incucyte image segmentation & analysis
%   - Segmentation and quantification of Incucyte myoblasts images 

%% Andra Chincisan, Institute of Neuropathology, USZ
% 2017

%% This code uses functions from Image Processing Toolbox MATLAB (%
% https://ch.mathworks.com/help/images/index.html ) and open-source code
% available on matlabcentral repository https://ch.mathworks.com/matlabcentral/

%% Initialize
clc
clear all
close all

%% Load Images
path = 'path\'
image_type = '.tif';
srcFiles = dir(strcat(path, '*', image_type));
index = 1;
count_all_images = 0;

for i = 1 : length(srcFiles)
    disp ('start ... ')
    close all;    
    image = strcat(path,srcFiles(i).name);  
    
    if srcFiles(i).name(8:9) == 'A1' | srcFiles(i).name(8:9) == 'B1' | srcFiles(i).name(8:9) == 'C1'  | srcFiles(i).name(8:9) == 'D1'  
        parameters_um = zeros (1,5);
        I = im2double(imread(image));        
        LfinalBW = zeros (size(I,1),size(I,2));
    else
        I = im2double(imread(image));
        original = im2double(imread(image));
        %% Segmentation  
        % Segmentation from gray k clustering 
        BW = adaptth_postpr_function(I); %Kmeans_image        
        %% Morphological operations 
        % Dilation
        BW_dil = dilation_function(BW);  
        % Fix holes
        BW_hol = fill_holes_areaop_function(BW_dil); 
        % Erosion
        BW = erosion_function(BW_hol);
        % Median filtering 
        Lfinal = medfilt2(BW);
        % Remove small obj 
        LfinalBW = remove_obj(Lfinal); % LfinalBW = remove_obj(Lfinal);
        % Calculate        
        [mean_parameters_um, parameters_um] = calculate_parameters(LfinalBW);
        %figure, imshowpair(original, LfinalBW, 'montage'); title('Final Segmentation');  
    end
    image_name = string (srcFiles(count_all_images + 1).name);
    namet =  {char(image_name)};
    [s1 s2 ] = size(parameters_um);
    count_all_images = count_all_images + 1
    if s1 > 1 
        no_fib = s1;
        all_param(count_all_images, :) = [sprintf('%d', count_all_images) namet sprintf('%f', mean_parameters_um(1)) sprintf('%f', mean_parameters_um(2)) sprintf('%f', mean_parameters_um(3)) sprintf('%f', mean_parameters_um(4)) sprintf('%f',mean_parameters_um(5)) sprintf('%d', no_fib)];   
        params_vector = [namet sprintf('%f', mean_parameters_um(1)) sprintf('%f', mean_parameters_um(2)) sprintf('%f', mean_parameters_um(3)) sprintf('%f', mean_parameters_um(4)) sprintf('%f',mean_parameters_um(5)) sprintf('%d', no_fib)];   
    else 
        no_fib = 0
        all_param(count_all_images, :) = [sprintf('%d', count_all_images) namet  sprintf('%f', parameters_um(1)) sprintf('%f', parameters_um(2)) sprintf('%f', parameters_um(3)) sprintf('%f', parameters_um(4)) sprintf('%f', parameters_um(5)) sprintf('%d', no_fib)];   
         params_vector = [namet  sprintf('%f', parameters_um(1)) sprintf('%f', parameters_um(2)) sprintf('%f', parameters_um(3)) sprintf('%f', parameters_um(4)) sprintf('%f', parameters_um(5)) sprintf('%d', no_fib)];   
    end

    index = index + 1;
    disp (all_param)
    imwrite(LfinalBW, strcat(path, '\Segmentation\', srcFiles(i).name, '_BW.png'));
end 
xlswrite(strcat (path , '\IncuCyteimages_results'), all_param)

%% --------------------------Functions ------------------------------------
%% Segmentation from gray k clustering and post-processing
% https://ch.mathworks.com/help/images/ref/adaptthresh.html
function BW = adaptth_postpr_function(image) %Kmeans_image
    % Red chanel : image(:,:,1)
    % Green : image(:,:,2)
    [s1, s2] = size(image);
    intensity_channel1 = sum (sum (image(:,:,1)))/ (s1*s2);
    intensity_channel2 = sum (sum (image(:,:,2)))/ (s1*s2);
    intensity_channel3 = sum (sum (image(:,:,3)))/ (s1*s2);
    if (intensity_channel1 > intensity_channel2) && (intensity_channel1 > intensity_channel3)
        Image_seg = image(:,:,1); 
         disp('red')
    end
    if (intensity_channel2 > intensity_channel1) && (intensity_channel2 > intensity_channel3)
        Image_seg = image(:,:,2); 
        disp('green')
    end
    if (intensity_channel3 > intensity_channel1) && (intensity_channel3 > intensity_channel2)
        Image_seg = image(:,:,3); 
    end
    % Adjust image contrast
    Image_adj = imadjust(Image_seg,[0.15 0.4],[]); %[0.01 0.4]
    % Median filtering
    Image_fil = medfilt2(Image_adj);
    % Adaptive threshold segmentation,TA - threshold, BW - segmented image
    TA = adaptthresh(Image_fil, 0.4);
    BW = imbinarize(Image_fil,TA);    
end

 %% Dilation line element
 % Two dilation operations with different structuring elements will be performed
 % https://ch.mathworks.com/help/images/index.html
function BW_dil = dilation_function(BW)  
    % First structuring element 
    se1 = strel('line', 10, 135);      
    % Second structuring element 
    se2 = strel('line', 10, 45); 
    % Operations
    BW_dil = imdilate(BW,se1); 
    BW_dil = imdilate(BW_dil,se2);  
    %figure, imshow(BW_dil);  title('dilation');
end

%% Fix holes
% https://ch.mathworks.com/help/images/ref/imfill.html
function BW = fill_holes_areaop_function(BW_dil) 
    % Close holes
    BWfill = imfill(BW_dil,'holes');
    % Area open (Remove objects smaller than 5000 pixels)
    BW = bwareaopen(BWfill, 20);
    %figure, imshow(BW), title('dilation + fix holes'); %title('bwareaopen'); 
end

%% Erosion
% Four erosion operatiosn with different structuring elements will be performed
% https://ch.mathworks.com/help/images/index.html
function BW = erosion_function(BW_hol)      
    % First structuring element 
    se1 = strel('line', 3, 0);
    % First structuring element 
    se2 = strel('disk',1);
    % Erosion
    BW = imerode(BW_hol, se1);
    BW = imerode(BW, se2);
    BW = imerode(BW, se2);    
    BW = imerode(BW, se2);     
end

%% Calculate
% https://ch.mathworks.com/help/images/ref/bwlabel.html
% https://ch.mathworks.com/help/images/ref/regionprops.html
function [mean_parameters_um, parameters_um] = calculate_parameters(Lfinal)
    % Parameters:
        % - Area
        % - length
        % - Perimeter
        % - Breath
        % - Shape
    
    % Stats
    if sum(sum(Lfinal)) > 0
        stats = regionprops('table',Lfinal, 'Area','MajorAxisLength','Perimeter' ,'BoundingBox');
        % Calculater values
        perimeter = stats.Perimeter;
        area = stats.Area;
        %MinorAxis = stats.MinorAxisLength;
        breadth = (perimeter - sqrt(power(perimeter,2) - 16*area))/4;
        shape = (4*pi*area./power(perimeter,2));
        heights = stats.BoundingBox(:,4);
        % Create vector
        parameters = [area perimeter breadth heights shape]; 

        % Transform values from pixels to um 
        perimeter_um = stats.Perimeter * 0.61;
        area_um = stats.Area * 0.61;
        breath_um = 0.61 * (perimeter - sqrt(power(perimeter,2) - 16*area))/4;
        shape_um = (4*pi*area./power (perimeter,2));
        heights_um = 0.61 *stats.MajorAxisLength;
        % Create vector
        parameters_um = [area_um perimeter_um breath_um heights_um shape_um];
        %disp (parameters_um)
        mean_parameters_um = zeros (1,5); 
        mean_parameters_um = mean (parameters_um); 
    else
        parameters_um = zeros (1,5); 
        mean_parameters_um = zeros (1,5); 
    end
end

%% Remove small obj 
% https://ch.mathworks.com/help/images/ref/bwareaopen.html
function LfinalBW = remove_obj(Lfinal)
    LfinalBW = bwareaopen(Lfinal, 20);
end      
