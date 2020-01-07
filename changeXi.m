%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% changeXi
% Melissa Mai
% 
% Updates the adhesion on each bead depending on the motion:
%   Trailing arm extension:     [high high low]
%   Leading arm extension:      [high high low]
%   Trailing arm contraction:   [low high high]
%   Leading arm contraction:    [low high high]
%
% INPUTS
% thisW                     Current motion [trailing leading]
% xi                        Global adhesion parameter
% xibounds                  [High Low] scale factors for adhesion
%
%
% OUTPUT
% thisXi                    1x3 vector of the adhesion values for the cell
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function thisxi = changeXi(thisW, xi, xibounds)
    high = xibounds(1);
    low = xibounds(2);
    if thisW(1) > 0
        % Trailing arm extension
        thisxi = [high high low];
    elseif thisW(2) > 0
        % Leading arm extension
        thisxi = [high high low];
    elseif thisW(1) < 0
        % Trailing arm contraction
        thisxi = [low high high];
    else
        % Leading arm contraction
        thisxi = [low high high];
    end
    
    thisxi = thisxi * xi;
end