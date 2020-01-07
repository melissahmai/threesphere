%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% racm3
% Melissa Mai
% Rotate a 3-sphere swimmer/crawler about its center of mass, 3d
%
% INPUTS
% coords                    3x3 (row corresponding to point)
% theta                     polar angle (wrt z-axis)
% phi                       azimuthal angle (wrt to x-axis in xy-plane)
%
% 
% OUTPUT
% newcoords                 3x3 matrix of coordinates for rotated cell
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newcoords = racm3(coords, theta, phi)

    % Principal orientation vector (new)
    alpha = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
    
    % Center of mass
    com = mean(coords);
    r0 = coords - com;
    
    % Current orientation
    ohat = (r0(3,:)-r0(1,:))/norm(r0(3,:)-r0(1,:));
    
    % Rotate by projecting onto alpha
    r(1,:) = dot(r0(1,:),ohat)*alpha';
    r(2,:) = dot(r0(2,:),ohat)*alpha';
    r(3,:) = dot(r0(3,:),ohat)*alpha';
    
    % Reposition
    newcoords = r + coords(2,:);
    
end