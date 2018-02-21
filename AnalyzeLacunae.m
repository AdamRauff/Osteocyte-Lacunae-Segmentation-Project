% AnalyzeLacunae.m

% v1 Authored by Adam Rauff & Chelsea Heveran

% Part 1
% Calculate moments of inertia
% and axes of each lacuna.

% Part 2
% Eigen decomposition - principal moments


% call on Get voxel location script to acuire necessary preceding
% information
GetVoxelLocations;

%% Part 1 - Calculate moments of interia
% pre-allocate structure
fieldVal =cell(TotLacNum,1);
MomInt = struct('theta',fieldVal,'phi',fieldVal, 'I', fieldVal, 'V', fieldVal, 'D', fieldVal);

for i = 1:TotLacNum
    
    % pre-allocating variables to store moments of intertia, eigenvalues, 
    % eigenvectors, and degrees (phi and theta)
    MomInt(i).I = zeros(3,3);  
    MomInt(i).V = zeros(3,3);
    MomInt(i).D = zeros(3,3);
    
    % calculating moments of inertia matrix
    for j = 1 : NumLacVox(i)
        
        % Axes are shifted for each object. Each moment is calculated with a
        % coordinate system where the origin is the COM of the lacuna
        MomInt(i).I(1,1) = MomInt(i).I(1,1) + ((maskVoxelLoc(i).LacCoord(j,2)- maskVoxelLoc(i).COM(1,2))^2+...
        (maskVoxelLoc(i).LacCoord(j,3)- maskVoxelLoc(i).COM(1,3))^2)...
        + 1/12*(y_dim^2+z_dim^2);
   
        MomInt(i).I(2,2) = MomInt(i).I(2,2) + ((maskVoxelLoc(i).LacCoord(j,1)- maskVoxelLoc(i).COM(1,1))^2+...
        (maskVoxelLoc(i).LacCoord(j,3)- maskVoxelLoc(i).COM(1,3))^2)...
        + 1/12*(x_dim^2+z_dim^2);
   
        MomInt(i).I(3,3) = MomInt(i).I(3,3) + ((maskVoxelLoc(i).LacCoord(j,1) - maskVoxelLoc(i).COM(1,1))^2+...
        (maskVoxelLoc(i).LacCoord(j,2) -  maskVoxelLoc(i).COM(1,2))^2)...
        + 1/12*(x_dim^2+y_dim^2);
   
        MomInt(i).I(1,2) = MomInt(i).I(1,2) + (maskVoxelLoc(i).LacCoord(j,1)-maskVoxelLoc(i).COM(1,1))...
        * (maskVoxelLoc(i).LacCoord(j,2)-maskVoxelLoc(i).COM(1,2));
    
        MomInt(i).I(1,3) = MomInt(i).I(1,3)+(maskVoxelLoc(i).LacCoord(j,1)-maskVoxelLoc(i).COM(1,1))...
        *(maskVoxelLoc(i).LacCoord(j,3)-maskVoxelLoc(i).COM(1,3));
    
        MomInt(i).I(2,3) = MomInt(i).I(2,3)+(maskVoxelLoc(i).LacCoord(j,2)-maskVoxelLoc(i).COM(1,2))...
        *(maskVoxelLoc(i).LacCoord(j,3)-maskVoxelLoc(i).COM(1,3));
    
    end
end

%% Part 2 - Principal Moments
for i = 1 : TotLacNum
    
    % multpily moments of intertia by voxel volume
    MomInt(i).I = MomInt(i).I*(x_dim*y_dim*z_dim);
    
    % Make sure all moments outside the diagonal are negatives
    MomInt(i).I(1,2) = -MomInt(i).I(1,2);
    MomInt(i).I(1,3) = -MomInt(i).I(1,3);
    MomInt(i).I(2,3) = -MomInt(i).I(2,3);
    
    % This matrix is symmetric 
    MomInt(i).I(2,1) = MomInt(i).I(1,2);
    MomInt(i).I(3,1) = MomInt(i).I(1,3);
    MomInt(i).I(3,2) = MomInt(i).I(2,3);
    
    % decompose moments tensor into eigenvalus and eigen vectors
    [MomInt(i).V,MomInt(i).D] = eig(MomInt(i).I);
    
    % The Eigen value matrix is a diagonal matrix (hence D)
    % these values, the eigen values, correspond to the momnets of the
    % object (lacuna) about the principle axes
    
    % the Eigenvector matrix holds unit vectors the correspond to the
    % direction of the principle aces
    
    % the smallest eigenvalue, D(1,1) corresponds to the smallest moment,
    % and hence the longest (major) princple axis. This means that V(:,1)
    % is our major axis direction
end