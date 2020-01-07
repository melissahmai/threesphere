%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Smat
% Melissa Mai
%
% Calculates the mobility matrix for regularized Stokeslets with a 
% wall boundary condition using the method of images, as described in 
% Ainley et al (2008). Only for one pair. Use with build_vmat for full
% system.
% Wall is assumed to be the xy-plane (z = 0).
%
%
% INPUTS
% Xe                        Coordinates of point at which velocity will be 
%                               evaluated
% Xs                        Coordinates of the Stokeslet (force)
% eps                       Epsilon, blob cutoff parameter
% eta                       Viscosity
%
%
% OUTPUT
% S                         Mobility matrix for s -> e
%                               ie, v(s->e) = S*fs
%
% SEE ALSO
%   - build_vmat
%   - vfield_S
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = Smat(Xe, Xs, eps, eta)

    % Epsilon ^2
    e2 = eps^2;
    
    % Distance above surface
    h = Xs(3);
    Xsim = Xs - 2*h*[0 0 1];
    
    % Distance Squared
    Xo = Xe - Xs;
    X = Xe - Xsim;
    
    ro2 = sum(Xo.^2);
    r2 = sum(X.^2);
    
    % Regularized distance
    roeps = sqrt(ro2 + e2);
    reps2 = r2+e2;
    reps = sqrt(reps2);
    reps3 = reps2*reps;
    
    % Velocity regularization coefficients
    % From first blob, psi = 15*e2 / (8*pi*re^7)
    H2o = 1/(8*pi*roeps^3);
    H2 = 1/(8*pi*reps3);
    
    % dH1 = H1o - H1
    dH1 = ((1/roeps - 1/reps)/(8*pi) + e2*(H2o-H2))*[1 0 0; 0 1 0; 0 0 1];
    
    % H1'(xim)/xim
    H1r = -h*(1 + 3*e2/(reps2)).*[-1 0 0; 0 1 0; 0 0 1]/(4*pi*reps3);
    
    % H2'(xim)/xim
    H2r = -3/(8*pi*reps^5);
    
    % From second blob, psi = 3*e2 / (4*pi*re^5)
    D2 = -3/(4*pi*reps^5);
    D1 = 1/(4*pi*reps3) + e2*D2;
    
    % Other constants for convenience
    h2zr = 2*h*H2r*X(3);
    D2H2h2 = D2*h^2 - H2 - h2zr;
    
    % X component
    Sx = [-2*h*X(3)    0   X(1)*(h - X(3))];

    % Y Component
    Sy = [0 0 X(2)*(h-X(3))];

    % Z Component
    Sz = [h h (2*h - X(3))] .* X;

    % Full matrix
    S = (2*H2*[Sx; Sy; Sz] + ...
        (X')*D2H2h2.*[1 1 -1].*X + ...
        dH1 + H2o*(Xo')*Xo + ...
        D1*h^2*[1 0 0; 0 1 0; 0 0 -1] + ...
        H1r.*X(3) ...
        )/eta;
end