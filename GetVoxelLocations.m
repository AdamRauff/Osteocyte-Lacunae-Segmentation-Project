%% GetVoxelLocations.m
% Finds number of mask voxels and records their (x,y,z) locations.
% Locations are defined as the center of the voxel in units of microns.

% prepare workspace
clc; close all; clear all


%% flags
    
% print the number of objects throughout each step of the segmentation
PRINT_NUM_OBJ = true; 

% display a graph (made with isosurface) of the segmented lacunae with
% arrows modeling their principle axes
DISP_MOMENT_VECT = false;

% The following flag should not can only be true if DISP_MOMENT_VECT is
% true
% set to true to plot lines of closest lacunae
LAC_DIST_LINE = false;

% save 3D lacunar image
LAC_3D_IM = false;%*********** still needs to be added

% if DISP_MOMENT_VECT == false && LAC_DIST_LINE == true
%   Throw error. This should never happen
% end
%% Segment Lacunae function

% call on Segment Lacunae function in order to obtain a segmented (and
% labeled image)
[mask, CCobj, VoxDim, FolderName, MouseName, IMClass, thresh]= Segment_Lacunae(PRINT_NUM_OBJ);

% This flag has been utilized by this point. clean up workspace
clear PRINT_NUM_OBJ

if strcmpi(mask, 'nada')
    %--------------------------
    % user has clicked cancel
    % quit out of program
    %--------------------------
end

%% Prepare image metric variables
% calculate image dimensions
maskSize = size(mask);

% total number of voxels
x_size = maskSize(1);
y_size = maskSize(2);
z_size = maskSize(3);

%voxel dimension (As read from microscope)
x_dim = VoxDim(1);
y_dim = VoxDim(2);
z_dim = VoxDim(3);

% total number of voxels present in the image stack
numTotalVoxels = x_size*y_size*z_size;

% total image volume
ImStackVolume = numTotalVoxels*(y_dim*x_dim*z_dim); 

% Total number of lacunae present in image
LabMask = struct2cell(CCobj);
TotLacNum = LabMask{3,1};     

%Vector that keeps track of Lacunae ID (Lacunar Arrangment)
LacArr = 1:1:TotLacNum; 

% Calculate number of mask voxels and get number corresponding to mask
numMaskVoxels = nnz(mask); %nnz stands for number of nonzeros

%total volume of all lacunae: multiply number of lacunae voxels by the size
% a voxel (micron^3)
TotLacsVol = numMaskVoxels*(y_dim*x_dim*z_dim); 

%percentage of volume occupied by lacunae (%)
perLacVol = TotLacsVol/ImStackVolume; 

% Lacunar density: 1 lacuna per bone volume (1/micron^3) 
lacDensity = (numTotalVoxels*(y_dim*x_dim*z_dim))/TotLacNum; 

% Consider using regionprops or getting the number of voxels of each
% lacunae from labelmatrix / bwconncomp

% pre allocate
NumLacVox = zeros(TotLacNum,1);

%% Calculate COM and all voxels of each lacunae
% pre-allocate structure to store COM and lacunae Coordinates of each lacunae
fieldVal =cell(TotLacNum,1);

% field that stores the coordinates of the center of mass of each lacunae (microns)
% (x, y, z) format, and a field that stores the location of each voxel of
% each lacunae, and marks surface voxels
maskVoxelLoc = struct('COM',fieldVal,'LacCoord',fieldVal); 

% for each lacuna, find number of voxels, store (x,y,z) coordinates of every belonging voxel, and
% determine surface area voxels
for i = 1 : TotLacNum
    % obtain linear indices of ith Lacuna in mask
    tempInds = cell2mat(CCobj.PixelIdxList(i));
    
    %stores the number of voxels of each lacunae
    NumLacVox(i) = length(tempInds);
    
    % convert indices to subscripts
    [row, col, slic] = ind2sub(maskSize,tempInds);
    
    % pre-allocate matrix of lacunae (x,y,z) locations of each voxel
    maskVoxelLoc(i).LacCoord = zeros(length(tempInds), 4);
    
        % scroll through every voxel that is found for the lacuna
        for j=1:length(row)
        % WHY?????????
        % Please find instance of use of this information stored
        % (coordinates in microns) and comment here!!!!!!!
        % 1) Center of Mass (in microns) - AnalyzeLacunae
        % 2) Calculation of moments of intertia - AnalyzeLacunae
        
            x = row(j)*x_dim-0.5*x_dim; % voxel x location in microns, center of voxel
            y = col(j)*y_dim-0.5*y_dim; % voxel y location in microns, center of voxel
            z = slic(j)*z_dim-0.5*z_dim; % voxel z location in microns, center of voxel
            maskVoxelLoc(i).LacCoord(j,:) = [x y z 0];

%--------------------------------------------------------------------------
% IMPORTANT NOTE, image coordinates are ordinaily flipped from matrix
% subscripts! (not here),
% That is, coordinates are given in (x, y), where subsripts are given in (row, col).
% the row specified the y coordinate, while column specifies the x
% coordinate
%--------------------------------------------------------------------------
            % Determine if point lies on the surface of lacunae by
            % evaluating neighboring voxels. if neighbor voxel == 0, voxel
            % is marked as surface voxel
            if mask(row(j)-1,col(j),slic(j)) == 0 || mask(row(j)+1,col(j),slic(j)) == 0 || ... 
               mask(row(j),col(j)-1,slic(j)) == 0 || mask(row(j),col(j)+1,slic(j)) == 0 ||...
               mask(row(j),col(j),slic(j)-1) == 0 || mask(row(j),col(j),slic(j)+1) == 0  
               %if the 4th coulmn has a 1 then the point is on the surface 
                maskVoxelLoc(i).LacCoord(j,4) = 1; 
            end            
        end
        % Center of mass in units of microns
        maskVoxelLoc(i).COM = mean(maskVoxelLoc(i).LacCoord(:,1:3)); 
end

% stores the volume of each lacunae (micron^3)
LacVol = NumLacVox*(y_dim*x_dim*z_dim); % this is an array

% for debugging purposes, ensure LacVol and NumLacVox are both arrays of
% size: (TotLacNum, 1)

% measuring the average volume of lacunae
avgLacVol = mean(LacVol);

% average standard deviation
stdLacVol = std(LacVol);

% clear unused variables for a more clear workspace
clear x y z i j k r row col slic tempInds


