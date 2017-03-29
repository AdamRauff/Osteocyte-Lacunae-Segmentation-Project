
% call on Analyze Lacunae script to acquire necessary preceding
% information
AnalyzeLacunae;

% Until this point, mask has been a labeled matrix, as outputed from
% Segment_Lacunae()
% here it gets redefined to be a binary unsigned 8-bit image, with 0 for background,
% and 255 for objects
tempMask = zeros(x_size, y_size, z_size, 'uint8');
tempMask(mask~=0) = 255;

% temporarily erase the mask variable
clear mask

% assignm mask to be tempMask. This ensures mask is of class unsigned 8-bit 
mask = tempMask;

% erase the tempMask variable
clear tempMask

if DISP_MOMENT_VECT == true
    
    figure('Name','3D Rendering Moment directions'); hold on;
    fv = isosurface(mask,0.5);              % insert any value between 0 and 1 for the isovalue.

    % isosurface likes gradients and evaluates difference in values.
    % Since we have mask, a 3D matrix with object labled from 1 to TotLacNum,
    % we want to display all objects and choose an isovalue that is smaller
    % than 1, and bigger than 0

    % note the color of the faces is given in RGB, can be changed per user
    % preference. currently displaying some shade of green. This color
    % preference can be added as a flag, or some user option defined earlier
    patch(fv,'FaceColor',[61/255 163/255 47/255], 'EdgeColor','none'); hold on;


    % this is a shortcut for view(-37.5, 30) which defines the view point in
    % polar coordinates from 0 view(azimuith, elevation)
    view(3);

    % Reset axes to give a scale that allows to see the arrows
    axis([y_size-1.5*y_size y_size+0.5*y_size x_size-1.5*x_size x_size+0.5*x_size z_size-1.5*z_size z_size+0.5*z_size]);

    % lighiting
    camlight;
    lighting gouraud; % interpolates lighting

    % perspective
    camproj perspective; % parallel lines drawn from perspective. more closely correlates with human vision.
                         % displays parallel lines as converging

    

    % label axis (note X, and Y are switched here!)
    xlabel('Y axis','FontSize', 24);
    ylabel('X axis','FontSize', 24);
    zlabel('Z axis','FontSize', 24);

    hold on;
end

% pre-allocation
X = zeros(TotLacNum,1); 
Y = zeros(TotLacNum,1);
Z = zeros(TotLacNum,1);
U = zeros(TotLacNum,1);
V = zeros(TotLacNum,1);
W = zeros(TotLacNum,1);

% compose COM matrix, and direction of principle axis vectors
for i = 1:TotLacNum
    % Center of Mass is stored in units of microns.
    % Divinding by the voxel length units gets the position on the matrix as a subscript
    % This way the arrows are plotted from the COM of each object
    X(i) = maskVoxelLoc(i).COM(2)/x_dim; 
    Y(i) = maskVoxelLoc(i).COM(1)/y_dim; 
    Z(i) = maskVoxelLoc(i).COM(3)/z_dim;
    
    U(i) = MomInt(i).V(2,1);  
    V(i) = MomInt(i).V(1,1);
    
    W(i) = MomInt(i).V(3,1);
    
    % This is a debugging trial
    % an attemp to flip the direction of the 3rd direction (z), and then
    % reset the vector (U, V, W) to be a unit vector
%     if sign(MomInt(i).V(3,1)) == -1
%         W(i) = -1-MomInt(i).V(3,1);
%     else
%         W(i) = 1-MomInt(i).V(3,1);
%     end
%     
%     U(i) = U(i)/norm(U(i));
%     V(i) = V(i)/norm(V(i));
%     W(i) = W(i)/norm(W(i));
end

if DISP_MOMENT_VECT == true
    
    % quiver plots a vector int the plot using the parameters for point and
    % direction specified in the above for loop
    quiver3(X,Y,Z,U,V,W,'b','AutoScale','on','AutoScaleFactor',0.8,'LineWidth',1.5); hold on;

    % the second quiver graphs the negative vecotr in order for the vector to
    % stretch across each lacunae
    quiver3(X,Y,Z,(-1*U),(-1*V),(-1*W),'r','AutoScale','on','AutoScaleFactor',0.8,'LineWidth',1.5); 
    
    
    % set the aspect ratio of the graph to the voxels dimensions
    daspect([1/VoxDim(1) 1/VoxDim(2) 1/VoxDim(3)]);
%     axis equal % axis equal sets the data aspect ratio as 1 to 1 to 1

    % the daspect ratio function sets the lengths of x, y, z specified as
    % appearing equal. We are looking to make each voxel look like an
    % elongated cuboid, with square subsection in x, y. 
    % Thus we set the lengths of the inverses of each dimension equal to each other. see documentation
    % for daspect()
end

% clear unnecessary variables to clean up workspace 
clear LabMats LabMask
