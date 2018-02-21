% LacunarDistance.m

% v1 Authored by Adam Rauff & Chelsea Heveran

% Compute distance of closest lacunae COM for each
% lacunae (COM to COM)

% Part 1
%   - Evaluate COM to COM distances for each lacuna

% Part 2
%   - find closest lacuna for each lacuna

% part 3 (optional flag)
%   - plot closest lacuna line

% part 4
%   - write data to csv

% call on orientation.m to run all necessary preceding scripts
orientation;

%% Part 1 - evaluate COM distances
% pre-allocate field in a structure
for i = 1:TotLacNum
    lacDist(i).COMDist = zeros(TotLacNum,1);
end
% ------------------------------------------------------------------------
% Note this field stores an array of size(TotLacNum,1) for each lacunae.
% Each array will hold the distance from that lacunae. 
% The array of lacunae 1 would read 0 in the first entry, as it is 0 microns away from its own
% COM, and then start listing magnitudes of distance in microns.
% ------------------------------------------------------------------------

% for each lacunae, scroll through all lacunae (other then itself)
for i = 1 : TotLacNum
    for j = 1 : TotLacNum
        if i ~=j
            
        % calcualte the distance from COM to COM using the 3D distance
        % formula: dst = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
        lacDist(i).COMDist(j,1) = sqrt((maskVoxelLoc(i).COM(1,1)-maskVoxelLoc(j).COM(1,1))^2 + ...
                                       (maskVoxelLoc(i).COM(1,2)-maskVoxelLoc(j).COM(1,2))^2 + ...
                                       (maskVoxelLoc(i).COM(1,3)-maskVoxelLoc(j).COM(1,3))^2);
        end
    end
end

%% Part 2 - find closest lacunae
% pre- allocate
shrtCOMDist = zeros(TotLacNum,1);

% This separate field in the same structure stores the closest lacunae by
% looking for the minimum value in the COMDist array that isn't zero
for i = 1:TotLacNum
    lacDist(i).shrtDist = min(lacDist(i).COMDist(lacDist(i).COMDist~=0)); % the second structure is just to keep track of the center of mass
    shrtCOMDist(i,1) = lacDist(i).shrtDist; % store in vector for ease of data conversion to csv
end

%% Part 3 - plotting
% if the flag for plotting the closest lacunae line is set to true (found
% in GetVoxelLocations.m script) then execute the following block
if LAC_DIST_LINE == true
    
    % pre - allocate variable
    k = zeros(TotLacNum,1);

    X = [0,0];
    Y = [0,0];
    Z = [0,0];

    % loop through each COMDist array to mark which lacunae is closet to each
    % lacunae. This is only used for plotting the shortest distance line
    for i = 1:TotLacNum
        for j = 1:TotLacNum
            if lacDist(i).COMDist(j) == lacDist(i).shrtDist
                k(i) = j;
            end
        end
    end


    % plot the line of shortest distance between lacunae in dashed magenta
    for i = 1:TotLacNum
        X = [maskVoxelLoc(i).COM(1)/x_dim, maskVoxelLoc(k(i)).COM(1)/x_dim]; % each coordinate is divided by its voxel dimension to convet form micron to voxels
        Y = [maskVoxelLoc(i).COM(2)/y_dim, maskVoxelLoc(k(i)).COM(2)/y_dim];
        Z = [maskVoxelLoc(i).COM(3)/z_dim, maskVoxelLoc(k(i)).COM(3)/z_dim];

        plot3(Y,X,Z,'-.m');
        hold on
    end
    
end

% clear workspace
clear i j LAC_DIST_LINE fieldVal CCobj DISP_MOMENT_VECT

%% Part 4 - Write data to CSV file Write 

% flip array to be column vector
if size(LacArr,1) == 1
    LacArr = LacArr';
end

% Data(1,1) = 
Titles = {'Lacunar ID', 'Ellipsoidal Surface Area', 'Ellipsodial Anisotropy', ...
             'Ellipsodial Volume','Lacunar Surface Area','Lacunar Anisotropy', ...
             'Lacunar Volume', 'Theta','Phi', 'Span Theta', 'Closest COM [micron]', ...
             'Long Positive Radius','Long Negative Radius','Second Positive Radius',...
             'Second Negative Radius','Short Positive Radius','Short Negative Radius', ...
             'Ellipsoidal Long Radius','Ellipsoidal Second Radius','Ellipsoidal Short Radius'};

% Data = cell((TotLacNum+1),20);
% for i = 1:20
%     Data{1,i} = Titles{i};
% end

Data = [LacArr, eliSurfArea, eliAni, eliVol, lacSurfArea, ...
       lacAni, LacVol, ThetaVec, PhiVec, spanTheta, shrtCOMDist, ...
       lngPosRad, lngNegRad, sndPosRad, sndNegRad, shrtPosRad, shrtNegRad, ...
       eliLngRad, eliSndRad, elishrtRad];

% write Data to csv file
csvwrite(strcat(FolderName,'/', MouseName, ' ','data.csv'), Data);

% total volume in image
Data2(1,1) = ImStackVolume;

% total number of voxels present in the image stack
Data2(1,2) = numTotalVoxels;

% total number of lacunae
Data2(1,3) = TotLacNum;

%total volume of all lacunae
Data2(1,4) = TotLacsVol;

% size of each voxel
Data2(1,5) = (x_dim*y_dim*z_dim);

% class of image (i.e 8-bit)
Data2(1,6) = IMClass;

% threshold used. (this number is multiplied by the mean of each image
% slice to determine the threshold, so it is a somewhat dynamic threshold
Data2(1,7) = thresh;

% Image depth
Data2(1,8) = z_size;

% write Data to csv file (second time)
csvwrite(strcat(FolderName,'/', MouseName, ' ','data2.csv'), Data2);

