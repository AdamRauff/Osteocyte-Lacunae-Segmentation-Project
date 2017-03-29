% anisotropy.m
% this code calculates the lacunar and ellipsodial anistropies.
% for ellipsoid, this is currently defined as the ratio of the smallest to
% biggest eigenvalues.
% for Lacunar, this is the ratio of the smallest to biggest radius

% call on graph lacunae in order to invoke all preceding functions
graphlacunae;

% pre-allocate structure
fieldVal =cell(TotLacNum,1);
I = struct('princI',fieldVal);

% organize principle moments of intertia to (3,1) vectors in field princI
for i = 1: TotLacNum
    I(i).princI = zeros(3,1);
    for j = 1:3         
        I(i).princI(j,1) = MomInt(i).D(j,j);
    end
end

% Solve for ellipsoid with identical moments of intertia
for i = 1:TotLacNum
    % Check Equations again, I believe it arrived from Mcreadie 2004, but I
    % also recall some alteration
    elipMat = [0 1 1; 1 0 1; 1 1 0];
    elipMat = (1/5)*LacVol(i)*elipMat;
    I(i).radii = elipMat\I(i).princI; % solving for the 3 by 1 vecotr [a;b;c] of radii 
    I(i).radii = sqrt(I(i).radii);    % of each lacunae that is estimated as an ellipsoid
end                                   

% pre- allocation
eliAni = zeros(TotLacNum,1);

% ellipsoidal anisotropy (anisotropy of the solved ellipsod)
for i = 1:TotLacNum
    eliAni(i) = min(I(i).radii)/max(I(i).radii);
end

% might want to consider goodness of fit of each ellipsoid

% Measuring Lacunar Anisotropy. 
% This measurement is taken by beginning at the COM, and proceeding in the
% direction of the principle axes, until a 0 is encountered. This way you
% measure the "Lacunar" Radius.

%-------------------------------------------------------------------------
% It is highly inefficient to have 6 separate for / while loops that do
% such similar things. consider consolidating into one Giant for/ while
% loop that calculates all six directions (3 axes, -/+) at once.
%------------------------------------------------------------------------

% measuring real lacunar radius about positive long axis 

% the following piece of code assumes that mask is a logical matrix that
% contains 1s and 0s
X1=ones(TotLacNum,1);  % location in matrix (mask) along the slope of eigenvector
Y1=ones(TotLacNum,1);
Z1=ones(TotLacNum,1);
oldX1=ones(TotLacNum,1);
oldY1=ones(TotLacNum,1);
oldZ1=ones(TotLacNum,1);

micX1 = zeros(1,3);     
micY1 = zeros(1,3);     
micZ1 = zeros(1,3);      
micX = zeros(1,3);              
micY = zeros(1,3);             
micZ = zeros(1,3); 

for i = 1:TotLacNum
t = 1;
while mask(round(Y1(i)),round(X1(i)),round(Z1(i))) ==  mask(round(oldY1(i)),round(oldX1(i)),round(oldZ1(i)))                                           
    % --------------------------------------------------------------------
    % the direction is measured from the moments that were calculated from
    % the micron values, shouldn't there be some conversion from the micron
    % to voxel?
    % --------------------------------------------------------------------
    X1(i) = X(i)+t*MomInt(i).V(2,1);           
    Y1(i) = Y(i)+t*MomInt(i).V(1,1);              
    Z1(i) = Z(i)+t*MomInt(i).V(3,1);
    
    oldX1(i) = X(i)+(t-1)*MomInt(i).V(2,1);
    oldY1(i) = Y(i)+(t-1)*MomInt(i).V(1,1);
    oldZ1(i) = Z(i)+(t-1)*MomInt(i).V(3,1);
    
    if round(X1(i)) == 0 || round(X1(i)) < 0
        X1(i) = 1; 
        oldX1(i) = 2;
        break
    elseif round(Y1(i)) == 0 || round(Y1(i)) < 0    % In order to avoid running into the extremities of the image
        Y1(i) = 1;
        oldY1(i) = 2;
        break
    elseif round(Z1(i)) == 0 || round(Z1(i)) < 0
        Z1(i) = 1;
        oldZ1(i) = 2;
        break
    elseif round(X1(i)) == x_size || round(X1(i)) > x_size
        X1(i) = x_size;
        oldX1(i) = x_size-1;
        break
    elseif round(Y1(i)) == y_size || round(Y1(i)) > y_size
        Y1(i) = y_size;
        oldY1(i) = y_size-1;
        break
    elseif round(Z1(i)) == z_size || round(Z1(i)) > z_size
        Z1(i) = z_size;
        oldZ1(i) = z_size-1;
        break
    end
    
    t=t+1;                                                                                                      
end


% get location of furthermost voxel and then measure distance using distance formula
micX1(i) = round(oldX1(i))*x_dim-0.5*x_dim;     % Not sure which point to use for the distance formula.
micY1(i) = round(oldY1(i))*y_dim-0.5*y_dim;     % should I use round() numbers since that marks the voxel location in the matrix?
micZ1(i) = round(oldZ1(i))*z_dim-0.5*z_dim;     % should I use the first location outside lacunae, last voxel in lacunae, 
micX(i) = X(i)*x_dim;              % or last voxel in lacunae with the center translated to middle of voxel?
micY(i) = Y(i)*y_dim;              % Right now it calculates the distance from COM of lacunae to location of first voxel outside lacunae
micZ(i) = Z(i)*z_dim; 
Rad(i).lngPos = sqrt((micX(i)-micX1(i))^2+(micY(i)-micY1(i))^2+(micZ(i)-micZ1(i))^2);
end

%measuring real lacunar radius about negative long axis

X1=ones(TotLacNum,1);  % location in matrix (mask) along the slope of eigenvector
Y1=ones(TotLacNum,1);
Z1=ones(TotLacNum,1);
oldX1=ones(TotLacNum,1);
oldY1=ones(TotLacNum,1);
oldZ1=ones(TotLacNum,1);

for i = 1:TotLacNum
t = 1;
while mask(round(Y1(i)),round(X1(i)),round(Z1(i))) ==  mask(round(oldY1(i)),round(oldX1(i)),round(oldZ1(i)))                                           

    X1(i) = X(i)-t*MomInt(i).V(2,1);           % direction with respect to the negative x axis
    Y1(i) = Y(i)-t*MomInt(i).V(1,1);              
    Z1(i) = Z(i)-t*MomInt(i).V(3,1);
    
    oldX1(i) = X(i)-(t-1)*MomInt(i).V(2,1);
    oldY1(i) = Y(i)-(t-1)*MomInt(i).V(1,1);
    oldZ1(i) = Z(i)-(t-1)*MomInt(i).V(3,1);
    
     if round(X1(i)) == 0 || round(X1(i)) < 0
        X1(i) = 1; 
        oldX1(i) = 2;
        break
    elseif round(Y1(i)) == 0 || round(Y1(i)) < 0    % In order to avoid running into the extremities of the image
        Y1(i) = 1;
        oldY1(i) = 2;
        break
    elseif round(Z1(i)) == 0 || round(Z1(i)) < 0
        Z1(i) = 1;
        oldZ1(i) = 2;
        break
    elseif round(X1(i)) == x_size || round(X1(i)) > x_size
        X1(i) = x_size;
        oldX1(i) = x_size-1;
        break
    elseif round(Y1(i)) == y_size || round(Y1(i)) > y_size
        Y1(i) = y_size;
        oldY1(i) = y_size-1;
        break
    elseif round(Z1(i)) == z_size || round(Z1(i)) > z_size
        Z1(i) = z_size;
        oldZ1(i) = z_size-1;
        break
    end
    
    t=t+1;                                                                                                      
end

% get location of furthermost voxel and then measure distance using distance formula
micX1(i) = round(oldX1(i))*x_dim-0.5*x_dim;     % Not sure which point to use for the distance formula.
micY1(i) = round(oldY1(i))*y_dim-0.5*y_dim;     % should I use round() numbers since that marks the voxel location in the matrix?
micZ1(i) = round(oldZ1(i))*z_dim-0.5*z_dim;     % should I use the first location outside lacunae, last voxel in lacunae, 
micX(i) = X(i)*x_dim;              % or last voxel in lacunae with the center translated to middle of voxel?
micY(i) = Y(i)*y_dim;              % Right now it calculates the distance from COM of lacunae to location of first voxel outside lacunae
micZ(i) = Z(i)*z_dim; 
Rad(i).lngNeg = sqrt((micX(i)-micX1(i))^2+(micY(i)-micY1(i))^2+(micZ(i)-micZ1(i))^2);
end

%measuring real lacunar radius about the second lonest positive axis

X1=ones(TotLacNum,1);  % location in matrix (mask) along the slope of eigenvector
Y1=ones(TotLacNum,1);
Z1=ones(TotLacNum,1);
oldX1=ones(TotLacNum,1);
oldY1=ones(TotLacNum,1);
oldZ1=ones(TotLacNum,1);

for i = 1:TotLacNum
t = 1;
while mask(round(Y1(i)),round(X1(i)),round(Z1(i))) ==  mask(round(oldY1(i)),round(oldX1(i)),round(oldZ1(i)))                                           

    X1(i) = X(i)+t*MomInt(i).V(2,2);           % direction with respect to the negative x axis
    Y1(i) = Y(i)+t*MomInt(i).V(1,2);              
    Z1(i) = Z(i)+t*MomInt(i).V(3,2);
    
    oldX1(i) = X(i)+(t-1)*MomInt(i).V(2,2);
    oldY1(i) = Y(i)+(t-1)*MomInt(i).V(1,2);
    oldZ1(i) = Z(i)+(t-1)*MomInt(i).V(3,2);
    
     if round(X1(i)) == 0 || round(X1(i)) < 0
        X1(i) = 1; 
        oldX1(i) = 2;
        break
    elseif round(Y1(i)) == 0 || round(Y1(i)) < 0    % In order to avoid running into the extremities of the image
        Y1(i) = 1;
        oldY1(i) = 2;
        break
    elseif round(Z1(i)) == 0 || round(Z1(i)) < 0
        Z1(i) = 1;
        oldZ1(i) = 2;
        break
    elseif round(X1(i)) == x_size || round(X1(i)) > x_size
        X1(i) = x_size;
        oldX1(i) = x_size-1;
        break
    elseif round(Y1(i)) == y_size || round(Y1(i)) > y_size
        Y1(i) = y_size;
        oldY1(i) = y_size-1;
        break
    elseif round(Z1(i)) == z_size || round(Z1(i)) > z_size
        Z1(i) = z_size;
        oldZ1(i) = z_size-1;
        break
    end
    
    t=t+1;                                                                                                      
end

% get location of furthermost voxel and then measure distance using distance formula
micX1(i) = round(oldX1(i))*x_dim-0.5*x_dim;     % Not sure which point to use for the distance formula.
micY1(i) = round(oldY1(i))*y_dim-0.5*y_dim;     % should I use round() numbers since that marks the voxel location in the matrix?
micZ1(i) = round(oldZ1(i))*z_dim-0.5*z_dim;     % should I use the first location outside lacunae, last voxel in lacunae, 
micX(i) = X(i)*x_dim;              % or last voxel in lacunae with the center translated to middle of voxel?
micY(i) = Y(i)*y_dim;              % Right now it calculates the distance from COM of lacunae to location of first voxel outside lacunae
micZ(i) = Z(i)*z_dim; 
Rad(i).sndPos = sqrt((micX(i)-micX1(i))^2+(micY(i)-micY1(i))^2+(micZ(i)-micZ1(i))^2);
end

%measuring real lacunar radius about the second lonest negative axis

X1=ones(TotLacNum,1);  % location in matrix (mask) along the slope of eigenvector
Y1=ones(TotLacNum,1);
Z1=ones(TotLacNum,1);
oldX1=ones(TotLacNum,1);
oldY1=ones(TotLacNum,1);
oldZ1=ones(TotLacNum,1);

for i = 1:TotLacNum
t = 1;
while mask(round(Y1(i)),round(X1(i)),round(Z1(i))) ==  mask(round(oldY1(i)),round(oldX1(i)),round(oldZ1(i)))                                           

    X1(i) = X(i)-t*MomInt(i).V(2,2);           % direction with respect to the negative x axis
    Y1(i) = Y(i)-t*MomInt(i).V(1,2);              %lacunae number 1, longest principal axis is on direction of first column of eigenvector
    Z1(i) = Z(i)-t*MomInt(i).V(3,2);
    
    oldX1(i) = X(i)-(t-1)*MomInt(i).V(2,2);
    oldY1(i) = Y(i)-(t-1)*MomInt(i).V(1,2);
    oldZ1(i) = Z(i)-(t-1)*MomInt(i).V(3,2);
    
    if round(X1(i)) == 0 || round(X1(i)) < 0
        X1(i) = 1; 
        oldX1(i) = 2;
        break
    elseif round(Y1(i)) == 0 || round(Y1(i)) < 0    % In order to avoid running into the extremities of the image
        Y1(i) = 1;
        oldY1(i) = 2;
        break
    elseif round(Z1(i)) == 0 || round(Z1(i)) < 0
        Z1(i) = 1;
        oldZ1(i) = 2;
        break
    elseif round(X1(i)) == x_size || round(X1(i)) > x_size
        X1(i) = x_size;
        oldX1(i) = x_size-1;
        break
    elseif round(Y1(i)) == y_size || round(Y1(i)) > y_size
        Y1(i) = y_size;
        oldY1(i) = y_size-1;
        break
    elseif round(Z1(i)) == z_size || round(Z1(i)) > z_size
        Z1(i) = z_size;
        oldZ1(i) = z_size-1;
        break
    end
    
    t=t+1;                                                                                                      
end

% get location of furthermost voxel and then measure distance using distance formula
micX1(i) = round(oldX1(i))*x_dim-0.5*x_dim;     % Not sure which point to use for the distance formula.
micY1(i) = round(oldY1(i))*y_dim-0.5*y_dim;     % should I use round() numbers since that marks the voxel location in the matrix?
micZ1(i) = round(oldZ1(i))*z_dim-0.5*z_dim;     % should I use the first location outside lacunae, last voxel in lacunae, 
micX(i) = X(i)*x_dim;              % or last voxel in lacunae with the center translated to middle of voxel?
micY(i) = Y(i)*y_dim;              % Right now it calculates the distance from COM of lacunae to location of first voxel outside lacunae
micZ(i) = Z(i)*z_dim; 
Rad(i).sndNeg = sqrt((micX(i)-micX1(i))^2+(micY(i)-micY1(i))^2+(micZ(i)-micZ1(i))^2);
end

%measuring real lacunar radius about the third lonest (shortest) positive axis

X1=ones(TotLacNum,1);  % location in matrix (mask) along the slope of eigenvector
Y1=ones(TotLacNum,1);
Z1=ones(TotLacNum,1);
oldX1=ones(TotLacNum,1);
oldY1=ones(TotLacNum,1);
oldZ1=ones(TotLacNum,1);

for i = 1:TotLacNum
t = 1;
while mask(round(Y1(i)),round(X1(i)),round(Z1(i))) ==  mask(round(oldY1(i)),round(oldX1(i)),round(oldZ1(i)))                                           

    X1(i) = X(i)+t*MomInt(i).V(2,3);           % direction with respect to the negative x axis
    Y1(i) = Y(i)+t*MomInt(i).V(1,3);              %lacunae number 1, longest principal axis is on direction of first column of eigenvector
    Z1(i) = Z(i)+t*MomInt(i).V(3,3);
    
    oldX1(i) = X(i)+(t-1)*MomInt(i).V(2,3);
    oldY1(i) = Y(i)+(t-1)*MomInt(i).V(1,3);
    oldZ1(i) = Z(i)+(t-1)*MomInt(i).V(3,3);
    
    if round(X1(i)) == 0 || round(X1(i)) < 0
        X1(i) = 1; 
        oldX1(i) = 2;
        break
    elseif round(Y1(i)) == 0 || round(Y1(i)) < 0    % In order to avoid running into the extremities of the image
        Y1(i) = 1;
        oldY1(i) = 2;
        break
    elseif round(Z1(i)) == 0 || round(Z1(i)) < 0
        Z1(i) = 1;
        oldZ1(i) = 2;
        break
    elseif round(X1(i)) == x_size || round(X1(i)) > x_size
        X1(i) = x_size;
        oldX1(i) = x_size-1;
        break
    elseif round(Y1(i)) == y_size || round(Y1(i)) > y_size
        Y1(i) = y_size;
        oldY1(i) = y_size-1;
        break
    elseif round(Z1(i)) == z_size || round(Z1(i)) > z_size
        Z1(i) = z_size;
        oldZ1(i) = z_size-1;
        break
    end
    
    t=t+1;                                                                                                      
end

% get location of furthermost voxel and then measure distance using distance formula
micX1(i) = round(oldX1(i))*x_dim-0.5*x_dim;     % Not sure which point to use for the distance formula.
micY1(i) = round(oldY1(i))*y_dim-0.5*y_dim;     % should I use round() numbers since that marks the voxel location in the matrix?
micZ1(i) = round(oldZ1(i))*z_dim-0.5*z_dim;     % should I use the first location outside lacunae, last voxel in lacunae, micX(i) = X(i)*x_dim;              % or last voxel in lacunae with the center translated to middle of voxel?
micY(i) = Y(i)*y_dim;              % Right now it calculates the distance from COM of lacunae to location of first voxel outside lacunae
micZ(i) = Z(i)*z_dim; 
Rad(i).shrtPos = sqrt((micX(i)-micX1(i))^2+(micY(i)-micY1(i))^2+(micZ(i)-micZ1(i))^2);
end

%measuring real lacunar radius about the third lonest (shortest) negative axis

X1=ones(TotLacNum,1);  % location in matrix (mask) along the slope of eigenvector
Y1=ones(TotLacNum,1);
Z1=ones(TotLacNum,1);
oldX1=ones(TotLacNum,1);
oldY1=ones(TotLacNum,1);
oldZ1=ones(TotLacNum,1);

for i = 1:TotLacNum
t = 1;
while mask(round(Y1(i)),round(X1(i)),round(Z1(i))) ==  mask(round(oldY1(i)),round(oldX1(i)),round(oldZ1(i)))                                           

    X1(i) = X(i)-t*MomInt(i).V(2,3);           % direction with respect to the negative x axis
    Y1(i) = Y(i)-t*MomInt(i).V(1,3);              %lacunae number 1, longest principal axis is on direction of first column of eigenvector
    Z1(i) = Z(i)-t*MomInt(i).V(3,3);
    
    oldX1(i) = X(i)-(t-1)*MomInt(i).V(2,3);
    oldY1(i) = Y(i)-(t-1)*MomInt(i).V(1,3);
    oldZ1(i) = Z(i)-(t-1)*MomInt(i).V(3,3);
    
     if round(X1(i)) == 0 || round(X1(i)) < 0
        X1(i) = 1; 
        oldX1(i) = 2;
        break
    elseif round(Y1(i)) == 0 || round(Y1(i)) < 0    % In order to avoid running into the extremities of the image
        Y1(i) = 1;
        oldY1(i) = 2;
        break
    elseif round(Z1(i)) == 0 || round(Z1(i)) < 0
        Z1(i) = 1;
        oldZ1(i) = 2;
        break
    elseif round(X1(i)) == x_size || round(X1(i)) > x_size
        X1(i) = x_size;
        oldX1(i) = x_size-1;
        break
    elseif round(Y1(i)) == y_size || round(Y1(i)) > y_size
        Y1(i) = y_size;
        oldY1(i) = y_size-1;
        break
    elseif round(Z1(i)) == z_size || round(Z1(i)) > z_size
        Z1(i) = z_size;
        oldZ1(i) = z_size-1;
        break
    end
    
    t=t+1;                                                                                                      
end

% get location of furthermost voxel and then measure distance using distance formula
micX1(i) = round(oldX1(i))*x_dim-0.5*x_dim;     % Not sure which point to use for the distance formula.
micY1(i) = round(oldY1(i))*y_dim-0.5*y_dim;     % should I use round() numbers since that marks the voxel location in the matrix?
micZ1(i) = round(oldZ1(i))*z_dim-0.5*z_dim;     % should I use the first location outside lacunae, last voxel in lacunae, 
micX(i) = X(i)*x_dim;              % or last voxel in lacunae with the center translated to middle of voxel?
micY(i) = Y(i)*y_dim;              % Right now it calculates the distance from COM of lacunae to location of first voxel outside lacunae
micZ(i) = Z(i)*z_dim; 
Rad(i).shrtNeg = sqrt((micX(i)-micX1(i))^2+(micY(i)-micY1(i))^2+(micZ(i)-micZ1(i))^2);
end

% pre - allocate
% placing the radii into arrays, each array represents a lacunae
lacRads = zeros(TotLacNum,6);

lngPosRad = zeros(TotLacNum,1);
lngNegRad = zeros(TotLacNum,1);
sndPosRad = zeros(TotLacNum,1);
sndNegRad = zeros(TotLacNum,1);
shrtPosRad = zeros(TotLacNum,1);
shrtNegRad = zeros(TotLacNum,1);

for i = 1:TotLacNum
    lacRads(i,:) = [(Rad(i).lngPos) (Rad(i).lngNeg) (Rad(i).sndPos) (Rad(i).sndNeg) (Rad(i).shrtPos) (Rad(i).shrtNeg)];
    
    % also store in vector for writing to csv later
    lngPosRad(i,1) = Rad(i).lngPos;
    lngNegRad(i,1) = Rad(i).lngNeg;
    sndPosRad(i,1) = Rad(i).sndPos;
    sndNegRad(i,1) = Rad(i).sndNeg;
    shrtPosRad(i,1) = Rad(i).shrtPos;
    shrtNegRad(i,1) = Rad(i).shrtNeg;
end

%measuring the anisotropy
lacAni = zeros(TotLacNum,1);
for i = 1:TotLacNum
    lacAni(i,:) = min(lacRads(i,:))/max(lacRads(i,:));
end

clear oldX1 oldY1 oldZ1 i j micX micX1 micY micY1 micZ micZ1