% Orientation.m

% v1 Authored by Adam Rauff & Chelsea Heveran

% This script computes the orientation of the principal axes of each lacuna
% with respect to the global coordinate system

% Theta - difference between primary principal axis of lacuna and x - axis
%       (projected onto 2D space - x,y)

% Phi - Difference between x-y plane and Z-axis

% Part 1 
%   - compute phi and theta

% Part 2
%   - compute statistics
%   - compute span theta

% call on surface area and all preceding code
surfaceArea;


%% Part 1
% pre-allocate
PhiVec = zeros(TotLacNum,1);
ThetaVec = zeros(TotLacNum,1);

% the structure MomInt has already been pre - allocated with the
% appropriate fileds in the AnalyzeLacune.m script

for i = 1:TotLacNum
    %obtaining theta for major axis
    MomInt(i).theta = atan2d(V(i),U(i));  
    MomInt(i).theta = MomInt(i).theta-90; % adjust, because x,y axis are flipped when consdiering image dimensions (coordinates and subscripts)
    
    % the following if statements are meant to portray theta in its
    % smallest magnitude. Since the direction of the principle axis in the
    % x, y plane can be theta, or opposite of theta, the following
    % operations adjust theta to be portrayed consistently for all lacunae.
    
    % These operations do not change the overall direction of the principle
    % axis
    
    % if theta is 180 or -180, adjust to 0
    if round(MomInt(i).theta) == 180 || round(MomInt(i).theta) == -180
        MomInt(i).theta = 0;
    end
    
    % if theta is really small, adjust to 0
    if abs(MomInt(i).theta) < 0.001 
        MomInt(i).theta = 0;
    end
    
    % if 90 < theta < 180, adjust theta to be its negative counterpart
    if MomInt(i).theta >= 90 && MomInt(i).theta <= 180
        MomInt(i).theta = 180-MomInt(i).theta;
    end
    
    % if -90 > theta > -180, adjust theta to be its positive counterpart
    if MomInt(i).theta <= -90 && MomInt(i).theta >= -180
        MomInt(i).theta = 180+MomInt(i).theta;
    end
    
    
    % if theta < -180, add 180 to make it a smaller number in
    % magnitude
    if MomInt(i).theta < -180
        MomInt(i).theta = MomInt(i).theta+180;
    end
    
    % if theta > 180, subtract 180 to make it smaller in magnitude
    if MomInt(i).theta > 180
        MomInt(i).theta = MomInt(i).theta-180;
    end
    
    MomInt(i).theta = -MomInt(i).theta;
    
    % store thetas in a vector. makes it easier to write to csv later
    ThetaVec(i,1) = MomInt(i).theta;
    
    MomInt(i).phi = acosd(W(i));

    if MomInt(i).phi == 180
        MomInt(i).phi = 0;
    end
end

%% Part 2
% standard summation by for loop 
avgTheta = 0;
for i = 1:TotLacNum
    avgTheta = avgTheta + MomInt(i).theta;
end
% find average by dividing by total number of entries 
% In this case the total number of lacunae
avgTheta = avgTheta/TotLacNum;

% standard summation by for loop 
avgPhi = 0;
for i = 1:TotLacNum
    avgPhi = avgPhi + MomInt(i).phi;
end

% find average
avgPhi = avgPhi/TotLacNum;

% Since we only care about direction, we are interested in getting the
% minimum (degree) value from the positive and negative directions
for i = 1:TotLacNum
    if MomInt(i).phi - avgPhi > acosd(-W(i)) - avgPhi
        MomInt(i).phi = acosd(-W(i));
    end
    PhiVec(i) = MomInt(i).phi;
end

% pre- allocate
VGenDirect = zeros(3,1);

% summing up all unit vectors to form the average direction of lacunae
for i = 1:TotLacNum
   VGenDirect = VGenDirect + [U(i); V(i); W(i)]; 
end

% find average
VGenDirect = VGenDirect/TotLacNum;

% normalize
VGenDirect = VGenDirect./norm(VGenDirect);

% pre- allocate
spanTheta = zeros(TotLacNum,1);

% meausring the angle between vectors. Angle only exists on the plane
% that is spanned by the two vectors
for i = 1:TotLacNum
    spanTheta(i) = acosd((dot(VGenDirect,[U(i); V(i); W(i)]))/(norm(VGenDirect)*norm([U(i); V(i); W(i)])));
    spanTheta(i) = abs(spanTheta(i)); % this should always result in positive number
    if spanTheta(i) > 90 && spanTheta(i) < 180
        spanTheta(i) = -(180-spanTheta(i)); 
    elseif spanTheta(i) > 180 
        spanTheta(i) = spanTheta(i)-180;
    end
end
