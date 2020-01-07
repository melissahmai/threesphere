%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% build_vmat
% Melissa Mai
%
% Build the velocity matrix for a three-sphere swimmer under the influence
% of a planar surface's hydrodynamics. 
% Hydrodynamics adapted from Ainley et al (2008)
% Wall is assumed to be the xy-plane (z = 0)
%
%
% INPUTS
% coords                    Coordinates of all three spheres, 3x3xncell
%                               [x1 y1 z1;
%                                x2 y2 z2;
%                                x3 y3 z3]
% ncells                    Number of cells
% eta                       Viscosity
% eps                       Epsilon, blob cutoff parameter (r)
% Xi                        1x(9*ncell) vector of xi values 
%                               [xi1x xi1y xi1z ... xi3y xi3z]
%
% OUTPUT
% vmat                      Full (9*ncell)x(9*ncell) matrix,
%                               [S1->1   S2->1  S3->1
%                                S1->2   S2->2  S3->2
%                                S1->3   S2->3  S3->3]
% 
%                           = [(s1x s1y s1z s2x s2y s2z s3x s3y s3z)->1x
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->1y
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->1z
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->2x
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->2y
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->2z
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->3x
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->3y
%                              (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->3z],
%
%                               where Sm->n is the mobility submatrix 
%                               describing the effect of particle m on n, 
%                               and (smi)->nj is the individual component
%                               describing the effect of particle m's 
%                               i-directional force on particle n's 
%                               j-directional velocity. Sn->n is the self 
%                               mobility tensor.
% vcube                     Deconstructed mobility matrix,
%                               3x3x(3*ncell)x(3*ncell). vcube(:,:,n,m) is
%                               S(m->n).
%
%
% SEE ALSO
%   - Smat
%   - calc_vf
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vmat, vcube] = build_vmat(coords, ncells, eta, eps, Xi)
    
    vmat = zeros(ncells*3);
    % Iterate through each sphere and build the submatrix
    for n = 1:(3*ncells)
        for m = 1:(3*ncells)
            % Insert submatrix into full vmat
            vmat(((n-1)*3 + 1):(n*3),((m-1)*3 + 1):(m*3)) = ...
                Smat(coords(n,:), coords(m,:), eps, eta);
        end
    end
    
    % Apply adhesion to get modified mobility matrix
    vmat = (eye(9*ncells) + vmat.*Xi)\vmat;
    vcube = permute(reshape(vmat,3,3*ncells,3,3*ncells), [1 3 2 4]);
end