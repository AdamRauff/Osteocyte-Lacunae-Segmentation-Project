%% Osteocyte isolation algorithm

% v1 Authored by Adam Rauff & Chelsea Heveran

% The purpose of this program is to segment osteocyte lacunae. 
% The segmentation is written for bone specimen stained with basic fuchsin,
% embedded in poly(methyl) methyl methacrylate, and polished to a thickness of 100-200
% microns.

% The images were acquired with a confocal microscope, and stored as 16
% bit tiff sequences

% Part 1
% data allocation
%   - read/write files
% elementary image provessing
%   - blurring, binarizing

% Part 2
% Opening - erosion followed by dilation

% Part 3
% Label each Lacunae (labelmatrix)
% obtain elementary statistics (volume)

% part 4
% volume filtering
% lacunae must be between 100 - 2000 micron^3 (line 246). These values can
% be changed.

% part 5
% Edge filtering - remove lacunae that touch borders

% part 6
% Save segmented image

% Part 7
% UI filtering - quality control by user

% Limitations of this Code: 
% the image set must consist of images of the same class
% and size. (class - 16 bit, 8 bit ...
%            size (Height x Width x Depth.
% Code written for tif image sets, all other formats are
% excluded (very easy to alter see line 107).
% Code assumes images are chronologically
% ordered for 3D stack alignment in the folder
% chosen above. (i.e Img001, Img002, ...)


function [UserFilteredMat, BWCCobj, VoxDim, FolderName, MouseName, IMClass, myThresh] = Segment_Lacunae(PRNT_NUM_OBJ)

%% Part1 - Data handling
global temp
temp = 0;

% Display a UI to allow user to select folder of interest
folder = uigetdir('','Enter the folder that contains images, then select Open');

% if folder is not found
if folder == 0
    disp('Program cancelled');
    UserFilteredMat = 'Nada';
    BWCCobj = 'Nada';
    VoxDim = 'Nada';
    FolderName = 'Nada';
    MouseName = 'Nada';
    return 
end

[~, MouseName, ~] = fileparts(folder);

% Need to insert GUI for entering image info: voxel size, slice thickness,
% slice size
prompt = {'x-dim (microns)','y-dim (microns)','slice thickness (microns)'};
      
dlg_title = 'Image Info';
answer = inputdlg(prompt,dlg_title);
S = size(answer);
if S(1) == 0 || S(2) == 0
    disp('Program Cancelled');
    return 
else  
    x_dim = str2double(cell2mat(answer(1))); % voxel dimension (x)
    y_dim = str2double(cell2mat(answer(2))); % voxel dimension (y)
    z_dim = str2double(cell2mat(answer(3))); % voxel dimension (z)
end

VoxDim = [x_dim, y_dim, z_dim];
% Create a new folder within the selected folder to store the processed
% images. Make Sure you can run code twice, without overwriting previous
% results
A = dir(strcat(folder,'/','Processed Image*'));
S = size(A);
if S(1) == 0;
    FolderName = [folder,'/','Processed Images ',date];
    BinTiffStack = [folder,'/','Processed Images ', date, ' Binary Tiffs'];
    mkdir(FolderName);
    mkdir(BinTiffStack); % make directory to store the binarized tiff images
else
    Folnum = S(1)+1;
    FolderName = [folder,'/','Processed Images',num2str(Folnum),' ',date];
    BinTiffStack = [folder,'/','Processed Images',num2str(Folnum),' ',date, ' Binary Tiffs'];
    mkdir(FolderName);
    mkdir(BinTiffStack);
end

% Get a list of all files in the folder (return a structure).
ImageFiles = dir(folder);

% pre-allocation of variables
tempNumImages = length(ImageFiles);
ListOfImageNames = cell(tempNumImages,1); % Initialize cell to store list of strings (image names)
filenames2 = cell(tempNumImages,1); % intialize cell to store list of strings (segmented image names)

% threshold used
myThresh = 1.9; % this number will me multiplied by the gaussian mean of each image to determine the threshold. This threshold worked for our data, and may need to change with different datasets.

% Filter list of files
for Index = 1:length(ImageFiles)
    
    % if image below certain size (10000 bytes), skip to next file.
    % Probably not 16-bit tiff.
    if ImageFiles(Index).bytes < 10000
        continue
    end
    
    % Get the base filename and extension.
    baseFileName = ImageFiles(Index).name;
    [~, name, extension] = fileparts(baseFileName);
    
    % Keep only TIF image files.
    if strcmpi(extension,'.tif')
                      
        temp = temp + 1; % variable that keeps track of the slice number
        
        ListOfImageNames{temp,1} = baseFileName;
        IM = imread(strcat(folder,'/',baseFileName));
        
        % check number of channels
        numChan = size(IM);
        
        if length(numChan) > 3 
            % exit
            disp('image contains too many dimensions');
            disp('Program cancelled');
            
            UserFilteredMat = 'Nada';
            BWCCobj = 'Nada';
            VoxDim = 'Nada';
            FolderName = 'Nada';
            MouseName = 'Nada';
            return
            
        elseif length(numChan) == 3 && numChan(3) > 1
            % exit
            disp('image contains too many color channels');
            disp('please provide grayscale images, 8 or 16 bit');
            
            UserFilteredMat = 'Nada';
            BWCCobj = 'Nada';
            VoxDim = 'Nada';
            FolderName = 'Nada';
            MouseName = 'Nada';
            return
        end

        % once the first image is read, pre-allocate 3D image Matrices
        if temp == 1

            imSize = size(IM); % size of every image
            imClass = class(IM); % class of every image
            % determine the maximum intensity of image class (used for
            % thresholding later)
            if strcmp(imClass, 'uint16')
                maxInt = 2^16-1;
                IMClass = 16;
            elseif strcmp(imClass, 'uint8')
                maxInt = 2^8-1;
                IMClass = 8;
            end
            OrigMat = zeros(imSize(1), imSize(2), tempNumImages, imClass); % 3D matrix that holds original images 
            SegMat3D = zeros(imSize(1), imSize(2), tempNumImages, imClass); % store segmented (binary) image in 3D matrix
        end
        OrigMat(:,:,temp) = IM; % 3D matrix that holds original images 
        filenames2{temp,1} = strcat(FolderName,'/',name,'PP','.tif'); %PP stands for Post Processed
        [SegIM] = Funrun(IM, maxInt, myThresh); % process image
        SegMat3D(:,:,temp) = SegIM; % store segmented (binary) image in 3D matrix
    end  
end

% remove unused entries in pre-allocated variables
i = tempNumImages;
boolVar = false;

while boolVar == false 
    if isempty(filenames2{i,1})
        filenames2(i) = [];
        ListOfImageNames(i) = [];
        OrigMat(:,:,i) = [];
        SegMat3D(:,:,i) = [];
    else
        boolVar = true;
    end
    
    i = i-1;
end

% filling all the holes in objects
SegMat3DF = imfill(SegMat3D,'holes');

%% Part 2 - Opening

% construct structuring element (3D sphere of radius 5)
SE = strel('sphere',3);

% perform erosion to get rid of attached canaliculi, then a dilation to
% dilate objects back to around their original volumes
Mat3D3 = imopen(SegMat3DF,SE);

%% Part 3 - Lable Lacunae 

CC2 = bwconncomp(Mat3D3);

L2 = labelmatrix(CC2); % labels each object found in the matrix

stats2 = regionprops(CC2,'Area');

%% Part 4 - Volume filtering

objVols = zeros(CC2.NumObjects,1); % pre-allocation of matrix

ind = 1;

if PRNT_NUM_OBJ == true
    disp(['perliminary # of objects: ', num2str(CC2.NumObjects)]);
end

% calculating volume of each object (method 1)
for i = 1:CC2.NumObjects
    
    % recall region props return a field named area, however this field
    % hold volume (number of voxels) when input argument is 3D
    objVols(i) = stats2(i).Area*(x_dim*y_dim*z_dim); % # of voxel * volume of each voxel (micron^3)
    
    % mark objects to filter by volume
    if objVols(i) < 100 || objVols(i) > 2000 % if the object is has less then 100 or more than 2000 micron^3, then it is assumed it is not a lacunae
        rmList(ind) = i;
        ind = ind + 1;
    end
end

L3 = L2;

% remove marked objects
if exist('rmList','var') 
    for i = 1:length(rmList)
        L3(L2==rmList(i)) = 0;
    end
end

if PRNT_NUM_OBJ == true
    CC3 = bwconncomp(L3);
    disp(['# of objects after volume filter: ', num2str(CC3.NumObjects)]);
end

clear rmList 

%% Part 5 - Edge filter
% Remove objects touching edges of image (cannot analyze these objects correctly 
% if the entire morphology isn't included)

% find indices of all non-zero voxels
Inds(:,1) = find(L3); % find non-zero elements
B = size(L3);

% get (x,y,z) subscripts of each nonzero element
[rows(:,1), cols(:,1), sli(:,1)] = ind2sub(B,Inds); 

clear Inds
rmList = []; % intialize

% filter through non-zero elements only to select objects that touch edge
% of matrix
for i = 1:length(rows)
    if rows(i) == 1 || rows(i) == B(1)
        rmList = [rmList; [rows(i), cols(i), sli(i)]];
    elseif cols(i) == 1 || cols(i) == B(2)
        rmList = [rmList; [rows(i), cols(i), sli(i)]];
    elseif sli(i) == 1 || sli(i) == B(3)
        rmList = [rmList; [rows(i), cols(i), sli(i)]];
    end
end

% convert subscripts back into indices (easier when calling on many matrix
% values)
Inds(:,1) = sub2ind(B, rmList(:,1), rmList(:,2), rmList(:,3));

% attain voxel depth values of marked object
objNums(:,1) = L3(Inds);

% remove duplicates and sort
objNums = unique(sort(objNums));

% remove marked objects
for i = 1:length(objNums)
    L3(L3==objNums(i)) = 0;
end

if PRNT_NUM_OBJ == true
    CC4 = bwconncomp(L3);
    disp(['# of objects after edge filter: ', num2str(CC4.NumObjects)]);
end

%% Part 6 - save segmented image superimposed on original

% convert back to binary 
BIM = L3>0.1;

% pre-allocate structure to store RGB images of segmented files
fieldVal =cell(B(3),1);
SegFileMat = struct('SegFile',fieldVal);

% pre-allocate matrix to store translucent segmented images
TransFileMat = zeros(B(1), B(2), B(3));  

% loop over each slice of image set
for i = 1:B(3)
    % super impose images
    TransFileMat(:,:,i) = imfuse(OrigMat(:,:,i),BIM(:,:,i),'diff');
    SegFileMat(i).SegFile = imoverlay(OrigMat(:,:,i),BIM(:,:,i),[1 0 0]);
end

%% Part 7 - Quality Control Check (GUI)

% create structure with information to send to GUI
UIStruct(1).ImageLoc = folder; % folder of original images
UIStruct(2).ImageLoc = ListOfImageNames; % list of names of original images within folder
UIStruct(3).ImageLoc = filenames2; % full path + image name of overlayed images

UIStruct(1).ImageSet = L3; % 3D matrix containing labled, mostly segmented objects
UIStruct(2).ImageSet = OrigMat; % 3D matrix containg the original images
UIStruct(3).ImageSet = SegFileMat; % 3D matrix containg the imoverlay images
UIStruct(4).ImageSet = TransFileMat; % 3D matrix containg the imfuse images

UIStruct(1).Properties = [x_dim, y_dim, z_dim]; % size of each voxel in microns

% call on GUI and pass structure
UIFiltMat = Image_Process_GUIDE(UIStruct);

% run the bwconncomp function and label matrix function again in order to
% produce chronocligally labeled objects
BWCCobj = bwconncomp(UIFiltMat);
UserFilteredMat = labelmatrix(BWCCobj);

disp(['# of objects after user filter: ', num2str(BWCCobj.NumObjects)]);

% write tif stack of binarized image to folder
for i = 1:B(3)
    imwrite(UserFilteredMat(:,:,i), [BinTiffStack, '/', 'Slice ', num2str(i), '.tif']);
end

end 

function [Binarized_image]  = Funrun(image, maxInt, myThresh)

% apply gaussian filter to image
Gauss_image = imgaussfilt(image,0.65);

% calculate the mean intensity of the image (used as threshold parameter
imGMean = mean2(Gauss_image);

% simple thresholding binarization
Binarized_image = im2bw(Gauss_image,(imGMean*myThresh)/maxInt); % this code runs the threshold as 1.9*mean of each image
end
