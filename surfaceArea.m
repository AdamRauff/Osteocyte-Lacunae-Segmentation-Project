% surfaceArea.m

% v1 Authored by Adam Rauff & Chelsea Heveran

% Measures surface area (SA) of each segmented lacunae, and 
% best fit ellipsoid 

% Part 1
% SA of ellipsoid

% Part 2
% SA from point cloud (delaunay triangulation)

% call on anisotropy and all preceding scripts 
anisotropy;

%% Part 1 - ellipsoid SA
% surface area of the ellipsoid
eliSurfArea = zeros(TotLacNum,1);

p = 1.6075;

% this is an equation for SA of ellipsoid
for i = 1:TotLacNum
    eliSurfArea(i) = 4*pi*((((I(i).radii(1))*(I(i).radii(2)))^(p)+((I(i).radii(1))*(I(i).radii(3)))^(p)+((I(i).radii(2))*(I(i).radii(3)))^(p))/3)^(1/p);
end

% pre- allocate
eliVol = zeros(TotLacNum,1);
eliLngRad = zeros(TotLacNum,1);
eliSndRad = zeros(TotLacNum,1);
elishrtRad = zeros(TotLacNum,1);

% measuring ellipsodial volume
for i = 1:TotLacNum
    eliVol(i,1) = 4/3*pi*I(i).radii(1)*I(i).radii(2)*I(i).radii(3);
    
    % place radii into vectors. used later when writing to a csv
    eliLngRad(i) = I(i).radii(1);
    eliSndRad(i) = I(i).radii(2);
    elishrtRad(i) = I(i).radii(3);
end

%% Part 2 - SA from point cloud
%finding location of surface points of each lacunae
% and store their subscripts
for r = 1:TotLacNum
    index = 1;
    for i = 1: NumLacVox(r)
        if maskVoxelLoc(r).LacCoord(i,4) == 1
            maskVoxelLoc(r).SurfCoord(index,:) = [maskVoxelLoc(r).LacCoord(i,1) maskVoxelLoc(r).LacCoord(i,2) maskVoxelLoc(r).LacCoord(i,3)];
            index = index +1;
        end
    end
end


%surface area from point cloud (delaunay Triangulation algorithm)
for i = 1:TotLacNum
    tri(i).tri = delaunayTriangulation(maskVoxelLoc(i).SurfCoord(:,1),maskVoxelLoc(i).SurfCoord(:,2),maskVoxelLoc(i).SurfCoord(:,3));
    [fbtri(i).fbtri, fbpoints(i).fbpoints] = freeBoundary(tri(i).tri);
end

lacSurfArea = zeros(TotLacNum,1);

% calculating area of a trinagle in 3D
 for r = 1:TotLacNum
    for i = 1:length(fbtri(r).fbtri(:,1))
        % Acquire first point of triangle
        firstP = [fbpoints(r).fbpoints(fbtri(r).fbtri(i,1),:)];
        firstP = firstP';
        
        % second point
        secP = [fbpoints(r).fbpoints(fbtri(r).fbtri(i,2),:)];
        secP = secP';
        
        % third point
        thirP = [fbpoints(r).fbpoints(fbtri(r).fbtri(i,3),:)];
        thirP = thirP';
        
        % calculate two of the vetices
        a = secP-firstP; 
        b = thirP-firstP;
        
        % sum area of all triangles of each lacuna
        lacSurfArea(r) = lacSurfArea(r)+(1/2)*norm(cross(a,b));
    end
 end

 % clear unnecessary variables from workspace
 clear fbpoints fbtri firstP tri i r j t a b index firstP secP thirP