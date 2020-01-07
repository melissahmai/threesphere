%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% calc_vf
% Melissa Mai
% 
% Calculates the forces and velocities for individual beads in the
% three-sphere model using a set of defined constraints:
% 
% For swimming:
%   Force-free (all three directions)
%   Torque-free (beta- and gamma- components)
%   Deformation velocities
%   Rigid body
% 
% For crawling:
%   Force-free (x & y)
%   Torque-free (beta)
%   Zero z-velocities
%   Deformation velocities
%   Rigid body
%
% Then calculates velocities using V = M*F, with M the mobility matrix
% (calcuated with build_vmat).
%
%
% INPUTS
% coords                    Coordinates for the cells, 3x3xncell
%                               [x1 y1 z1;
%                                x2 y2 z2;
%                                x3 y3 z3]
% eps                       Epsilon (regularization scale) / bead radius
% eta                       Viscosity
% thisW                     1x2xncell vector of the arm velocities
%                               [trailing leading]
% thisXi                    ncellx3 matrix of xi values
% fthresh                   Threshold force
%
%
% OUTPUTS
% V                         3x3xncell matrix of velocity components
% F                         3x3xncell matrix of internal force components
% theta                     1xncell vector of polar angles
% phi                       1xncell vector of azimuthal angles
% scale                     Threshold scaling factor, 
%                               min(1, fthresh / max(|F|))
% Fadh                      3x3xncell matrix of adhesive force components
%
% 
% SEE ALSO
%   - build_vmat
%   - Smat
%   - sphere3
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, F, theta, phi, scale, Fadh] = calc_vf(coords, eps, eta, ...
            thisW, thisxi, fthresh)
    
    % Number of cells
    ncells = size(coords,3);
    
    % Transform into a (ccells*3)x3 matrix for easier manipulation in vmat
    % calculation
    coordsT = reshape(permute(coords, [1 3 2]), ncells*3, 3);
    
    % Xi vector 
    Xi = zeros(1,9*ncells);
    
    for i = 1:ncells
        Xi((((i-1)*9)+1):(((i-1)*9)+2)) = thisxi(i,1); 
        Xi((((i-1)*9)+4):(((i-1)*9)+5)) = thisxi(i,2); 
        Xi((((i-1)*9)+7):(((i-1)*9)+8)) = thisxi(i,3);
    end
    
    % Build mobility matrix
    [m, M] = build_vmat(coordsT, ncells, eta, eps, Xi);
    
    % Forcefree template
    forcefree = [1 0 0 1 0 0 1 0 0;
        0 1 0 0 1 0 0 1 0;
        0 0 1 0 0 1 0 0 1];
    
    % Initialize
    Cmat = zeros(9*ncells);
    Cvec = zeros(9*ncells,1);
    phi = zeros(1, ncells);
    theta = zeros(1, ncells);
    
    % Build constraints for each cell
    for n = 1:ncells
        % Arm deformation velocity constraints
        % Leading arm
        Ll = coords(3,:,n) - coords(2,:,n);
        ll = sqrt(sum(Ll.^2));
        Wl = thisW(1,2,n);
        
        % Trailing arm
        Lt = coords(1,:,n) - coords(2,:,n);
        lt = sqrt(sum(Lt.^2));
        Wt = thisW(1,1,n);
        
        % Swimming angle and relevant basis vectors
        theta(n) = acos(Ll(3)/ll);
        phi(n) = atan(Ll(2)/Ll(1));
        
        % Correct for backward-swimmers
        fwd = (coords(3,1,n) - coords(1,1,n) > 0);
        if ~fwd 
            phi(n) = phi(n) + pi;
        end
        
        % Orthonormal orientation basis
        alpha = [sin(theta(n))*cos(phi(n)); 
            sin(theta(n))*sin(phi(n));
            cos(theta(n))];
        beta = [cos(phi(n))*cos(theta(n)); 
            sin(phi(n))*cos(theta(n)); 
            -sin(theta(n))];
        gamma = [-sin(phi(n)); cos(phi(n)); 0];
        
        
        
        % Torque-free constraint
        torquefree = ...
            [0 -Lt(3) Lt(2) 0 0 0 0 -Ll(3) Ll(2);
            Lt(3) 0 -Lt(1) 0 0 0 Ll(3) 0 -Ll(1);
            -Lt(2) Lt(1) 0 0 0 0 -Ll(2) Ll(1) 0];
        
        trow = [sum(torquefree .* beta); sum(torquefree .* gamma)];
        
        % Initializing
        wlrow = zeros(1,9*ncells); wtrow = zeros(1,9*ncells);
        brow = zeros(1,9*ncells); grow = zeros(1,9*ncells);
        
        % Iterate through each bead to fill in the rows for distortion
        % velocity and rigid body constraints
        for i = 1:(3*ncells)
            % Leading arm
            lead = (M(:,:,(n-1)*3+3,i) - M(:,:,(n-1)*3+2,i));
            
            % Trailing arm
            trail = (M(:,:,(n-1)*3+2,i) - M(:,:,(n-1)*3+1,i));
            
            % Convenient definition of relevant indices
            inds = (3*(i-1)+1):(3*i);
            
            % Projections onto alpha (W constraint)
            % Leading arm forces
            wlrow(inds) = sum(lead.*alpha);
            
            % Trailing arm forces
            wtrow(inds) = sum(trail.*alpha);
            
            % Projections onto beta (rigid body)
            brow(inds) = sum((lt*lead - ll*trail).*beta);
            
            % Projections onto gamma (rigid body)
            grow(inds) = sum((lt*lead - ll*trail).*gamma);
        end
        
        % Combine all rows
        Cmati = zeros(9, 9*ncells);
        Cmati(1:3, ((n-1)*9+1):(n*9)) = forcefree;
        Cmati(4:5, ((n-1)*9+1):(n*9)) = trow;
        
        Cmati(6:end, :) = [wlrow; wtrow; brow; grow];
        
        % Constraint value vector
        Cveci = [0 0 0 0 0 Wl Wt 0 0]';
        
        
        % Include z-velocity constraint if attached to the wall
        if (any(thisxi(n,:) ~= 0))
            % Release z-force, gamma torque constraint, theta constraint
            % Add on z-velocity constraint
            Cmati = [Cmati([1:2 4 6:7 9],:); m(((n-1)*9+3):3:(n*9),:)];
            Cveci = [Cveci([1:2 4 6:7 9]); 0; 0; 0];
        end
        
        % Add submatrix to full constraint matrix; same for vector
        Cmat(((n-1)*9+1):(n*9),:) = Cmati;
        Cvec(((n-1)*9+1):(n*9)) = Cveci;
    end
    
    
    % Solve for forces
    F = linsolve(Cmat, Cvec);
    f = reshape(F, 3, 3*ncells)';
    
    % Scale forces to threshold
    scale = min(1, fthresh / max(sqrt(sum(f.^2,2))));
    F = F * scale;
    
    % Solve for velocities
    V = m*F;
    
    % Adhesive
    Fadh = -Xi' .* V;
    
    % Reformat
    F = permute(reshape(F,3,3,ncells), [2 1 3]);
    V = permute(reshape(V,3,3,ncells), [2 1 3]);
    Fadh = permute(reshape(Fadh,3,3,ncells), [2 1 3]);
end