%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% sphere3
% Melissa Mai
%
% Rigid three-sphere swimmer/crawler (no internal bending), with 
% zero-z-velocity constraint when xi is active. Three-sphere swimmer based 
% off Najafi & Golestanian (2004). Surface constraints based off 
% Daddi-Moussa-Ider et al (2018). Regularized Blake tensors derived in
% Ainley et al (2018).
%
% This model uses discrete, constant-velocity deformation phases. Flow
% field calculations are not included in this script.
%
%
% INPUTS
% 'L', 'arm', 'length'      Mean length of each arm. Extensions and
%                               contractions, in addition to initial
%                               configuration, will be determined to
%                               maintain a mean arm length of L.
%                               Default: 1
% 'dt'                      Initial time step
%                               Default: 1e-2
% 'r', 'a'                  Radius of the beads. Assumes the same radius
%                               for all beads.
%                               Default: 0.1*L
% 'eta', 'viscosity'        Viscosity
%                               Default: 1
% 'h', 'height', 'z0',      Initial height of the swimmer's center of mass
%   'heights', 'wall'           off the wall (which is assumed to be the 
%                               xy-plane)
%                               Default: r
% 'com', 'cm', 'init_com',  Initial center of mass. Can be given as a
%   'initcom', 'config',        2-vector to specify projection onto the
%   'initconfig',               xy-plane and use the 'height' parameter for
%   'init_config'               the z-component, or can be given as a
%                               3-vector to override the 'height'
%                               parameter.
%                               Default: [0 0 height]
% 'ncycle', 'ncycles'       Total number of swimming cycles
%                               Default: 10
% 'xi', 'adh', 'adhesion'   Global adhesive friction parameter. Will be
%   'friction'                  scaled at each phase according to the
%                               movement:
%                               Leading+    [high high low]
%                               Leading-    [low high high]
%                               Trailing+   [high high low]
%                               Trailing-   [low high high]
%                               Change 'xibounds', 'xihigh', and/or 'xilow'
%                               to alter these scales
%                               Default: 0
% 'u0', 'dv', 'w'           Deformation/distortion velocity for the arms.
%                               Will use the same value of u0 for each
%                               deformation; to change this, see 'seq'.
%                               Default: 0.5*L
% 'deform'                  Length of full contraction/extension. For
%                               example, a 'deform' of 0.5 with a mean arm
%                               length L of 10 leads to arm lengths ranging
%                               from 7.5 to 12.5.
%                               Default: 0.5*L
% 'seq', 'sequence'         Sequence of arm motions. Typically given in
%                               terms of +/- 1 but can be changed to
%                               modulate the deformation velocity of each
%                               step. Each row represents a phase, and the
%                               columns represent [trailing leading] arms.
%                               Use +1 for extension and -1 for
%                               contraction. The initial configuration will
%                               be determined so that the first step can be
%                               completed.
%                               Default: [1 0; # Trailing extension
%                                         0 1; # Leading extension
%                                        -1 0; # Trailing contraction
%                                         0 -1]# Leading contraction
% 'phase', 'startphase',    Phase offset, in the range [0 4] (though if
%   'phaseoffset',              given as a value greater than four, it will
%   'offset'                    be rescaled via modulo to be within that
%                               range. Can be given as a vector to give
%                               individual offset values to each cell or as
%                               a single value to assign a global phase
%                               offset. The offset will determine the
%                               starting configuration of the cell, with
%                               [0 1) corresponding to the first phase,
%                               [1 2) corresponding to the second phase,
%                               etc. For example, an offset of 0.5 will
%                               start the swimmer halfway through the first
%                               phase defined by seq (see above). An offset
%                               of 2 will begin the cell at the start of
%                               the third motion phase.
%                               Default: 0
% 'fthresh', 'thresh'       Threshold force. All forces will be scaled if
%                               any of the required forces (the norms, not
%                               the components) exceeds the threshold.
%                               Forces are scaled as 
%                                   F * min(1, fthresh/max(fnorms))
%                               Default: 100
% 'theta', 'polar'          Polar angle (angle wrt z-axis). 
%                               Default: pi/2
% 'phi', 'azimuthal',       Azimuthal angle (angle in xy-plane wrt to
%   'azimuth'                   x-axis)
%                               Default: 0
% 'collide', 'collision'    Collision threshold. If any sphere approaches
%                               within 'collide' of the wall, the
%                               simulation stops. A notification of the
%                               collision and the step at which it occurs
%                               will be generated in the command line.
%                               Default: 0.75*r
% 'xihigh', 'high', 'xih'   Upper bound for the adhesion scaling (see
%                               'xi'). Overridden if 'xibounds' is defined.
%                               Default: 1
% 'xilow', 'low', 'xil'     Lower bound for the adhesion scaling (see
%                               'xi'). Overridden if 'xibounds' is defined.
%                               Default: 0.2
% 'xibounds'                Upper and lower bounds for adhesion scaling
%                               (see 'xi'). Of the form [high low], though
%                               it will be rearranged if [low high] is
%                               given. Overrides 'xihigh' and 'xilow' if
%                               given, otherwise, will be defined as
%                               [xihigh xilow].
%                               Default: [xihigh xilow]
% 'fiber'                   Boolean expression whether the simulation is
%                               run on a fiber (ie, substrate hydrodynamics
%                               are turned off by setting cell height
%                               arbitrarily large (r*1e4), but adhesion
%                               still exists. Height will be reset back to 
%                               height at the end of the simulation.) If 
%                               'height' parameter is explicitly defined, 
%                               the cell will still be detached from the 
%                               wall, but at the specified height (and not 
%                               arbitrarily far away)
%                               Default: false
% 'recenter'                Recenter the starting configuration so it is
%                               centered at the xy origin
%                               Default: false;
% 'tol', 'tolerance'        Tolerance for boundary cutoff; the motion phase
%                               will advance when the arm length is within
%                               tol of the actual value
%                               Default: 1e-6
% 'rdt', 'readout'          Readout time step. The program will not record
%                               every step (for outputs like allconf, V, F,
%                               thetavec, phivec, psivec, Lvec, Fadh, tvec)
%                               but will report in intervals of (at least)
%                               rdt.
%                               Default: 100*dt
%
% NOTE: Input names are not case-sensitive.
% NOTE: Presets of inputs can be defined and passed as a single cell: see
%       example. Extra inputs can be added after the preset, but only one
%       preset can be defined.
%
%
% OUTPUTS
% econf                     (ncycle+1)x7xncell matrix of the time, 
%                               coordinates of the middle sphere, and 
%                               angles at the end of each cycle. If the 
%                               cell collides with the wall (see 'collide' 
%                               above), this matrix will be truncated to 
%                               only include completed cycles. The matrix 
%                               includes the initial position at t=0.
%                               [time x2 y2 z2 theta phi psi]
% allconf                   (3*nstep+3)x3xncell matrix of the coordinates 
%                               at each time step. Coordinates at each time
%                               step are each triplet of rows (each row is
%                               one sphere).
% V                         (3*nstep+3)x3xncell matrix of the velocity 
%                               components at each time step. Same
%                               dimensions as allconf.
% F                         (3*nstep+3)x3xncell matrix of the force 
%                               components at each time step. Same
%                               dimensions as allconf.
% thetavec                  (nstep+1)xncell vector of the polar angles at 
%                               each time step
% phivec                    (nstep+1)xncell vector of azimuthal angles
% psivec                    (nstep+1)xncell vector of internal angles
% Lvec                      (nstep+1)x2xncell matrix of arm lengths. 
%                               Columns given as [trailing leading]
% tvec                      (nstep+1)x1 vector of time points
% cycvec                    (nstep+1)x1 vector of time points transformed
%                               to their relative position within the cycle
% hitthresh                 Boolean of if the threshold force (fthresh) was
%                               hit
%
% 
% USAGE
% Ex: >> sphere3
%   Will run the simulation with all default parameters
%
% Ex: >> sphere3('L', 20, 'height', 4, 'phi', pi/2)
%   Will use a mean arm length of 20 and will start the swimmer 4 units
%   above the wall with azimuthal angle pi/2 (swimming in the +y direction)
%
% Ex: >> preset1 = {'initconfig', [0 0 0; 1 3 0], 'phi', [0 pi], ... 
%                   'ncycle', 200, 'dt', 1e-3, 'rdt', 1e-1}
%     >> sphere3(preset1)
%   Allows pre-defined sets of inputs to be passed as a single input. This
%   example initializes two anti-parallel cells, one positioned at [0 0 0]
%   moving in the +x direction and the other at [1 3 0] moving in the -x
%   direction. Runs for 200 cycles with dt = 1e-3 and readout time step of
%   1e-1. Additional inputs can be included after the preset; ie, 
%       >> sphere3(preset1, 'L', 10)
%   is permissible, but only one preset is allowed. 
%       >> sphere3(preset1, preset2)
%   is not allowed.
%
% 
% SEE ALSO
%   - calc_vf
%   - racm3
%   - changeXi
%   - build_vmat
%   - Smat
%   - vcurve
%   - phicurve
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econf, allconf, V, F, Fadh, thetavec, phivec, psivec, ...
            Lvec, tvec, cycvec, hitthresh] ... 
            = sphere3(varargin)
    %% Input processing
    % If the first varargin is a cell (ie, preset), unpack
    if iscell(varargin{1})
        varargin = horzcat(reshape(varargin{1},1,length(varargin{1})), ...
            varargin{2:end});
    end

    fiber = false;

    for i = 1:2:length(varargin)
        this = varargin{i+1};
        switch lower(varargin{i})
            case {'l', 'arm', 'length'}
                L = this;
            case 'dt'
                dt = this;
            case {'r', 'a'}
                r = this;
            case {'eps', 'epsilon', 'blob'}
                eps = this;
            case {'eta', 'viscosity'}
                eta = this;
            case {'h', 'z0', 'height', 'wall', 'heights'}
                height = this;
                if strcmpi(height, 'none')
                    height = 1e10;
                end
            case {'com', 'cm', 'init_com', 'initcom', ...
                    'initconfig', 'config', 'init_config'}
                com = this;
            case {'ncycle', 'ncycles'}
                ncycle = this;
            case {'xi', 'adh', 'adhesion', 'fric'}
                xi = this;
            case {'u0', 'dv', 'w'}
                u0 = this;
            case {'deform'}
                deform = this;
            case {'seq', 'sequence'}
                seq = this;
            case {'fthresh', 'thresh'}
                fthresh = this;
            case {'theta', 'polar'}
                theta = this;
            case {'phi', 'azimuth', 'azimuthal'}
                phi = this;
            case {'collide', 'collision'}
                collide = this;
            case {'xihigh', 'high', 'xih'}
                xihigh = this;
            case {'xilow', 'low', 'xil'}
                xilow = this;
            case 'xibounds'
                xibounds = this;
            case 'fiber'
                fiber = this;
            case 'recenter'           
                recenter = this;
            case {'phase', 'startphase', 'phaseoffset', 'offset'}
                offset = this;
            case {'tol', 'tolerance'}
                tol = this;
            case {'rdt', 'readout'}
                Rdt = this;
        end
    end

    % Default values
    if ~exist('L', 'var'); L = 1; end
    if ~exist('dt', 'var'); dt = 1e-2; end
    if ~exist('r', 'var'); r = 0.1*L; end
    if ~exist('eps', 'var'); eps = r; end
    if ~exist('eta', 'var'); eta = 1; end
    if ~exist('height', 'var'); height = r; heightspec = false; 
    else; heightspec = true; end
    if ~exist('com', 'var'); com = [0 0 height]; end
    if ~exist('ncycle', 'var'); ncycle = 10; end
    if ~exist('xi', 'var'); xi = 0; end
    if ~exist('u0', 'var'); u0 = 0.1*L; end
    if ~exist('deform', 'var'); deform = 0.5*L; end
    if ~exist('seq', 'var'); seq = [1 0; 0 1; -1 0; 0 -1]; end
    if ~exist('fthresh', 'var'); fthresh = 1e5; end
    if ~exist('theta', 'var'); theta = pi/2; end
    if ~exist('phi', 'var'); phi = 0; end
    if ~exist('collide', 'var'); collide = 0.75*r; end
    if ~exist('xihigh', 'var'); xihigh = 1; end
    if ~exist('xilow', 'var'); xilow = 0.2; end
    if ~exist('xibounds', 'var'); xibounds = [xihigh xilow]; end
    if ~exist('recenter', 'var'); recenter = false; end
    if ~exist('offset', 'var'); offset = 0; end
    if ~exist('Rdt', 'var'); Rdt = 100*dt; end
    if ~exist('tol', 'var'); tol = 1e-6; end

    %% Useful parameters
    % Number of cells
    ncells = size(com,1);

    % Xi bounds
    xibounds = sort(xibounds, 'descend');

    % Deformation velocity sequence [trail lead]
    % Adjust sequence to number of cells
    if size(seq,3) ~= ncells
        seq = ones(4,2,ncells).*seq;
    end

    % Same with u0
    if size(u0,3) ~= ncells
        u0 = ones(1,1,ncells).*u0;
    end
    
    % Same with offset
    if size(offset) ~= ncells
        offset = ones(ncells,1).*offset(1);
    end

    % Deformation velocities
    W = seq.*u0;

    % Time and number of steps for each phase
    phasetime = deform/u0(1);

    % Time and number of steps for each cycle
    cyctime = phasetime * 4;


    %% Define initial coordinates for the swimmer
    if size(com, 2) == 2
        com = [com height];
    end
    % Initialize
    coords = zeros(3,3,ncells);
    h0 = zeros(ncells,1); g0 = zeros(ncells,1);

    if length(phi) ~= ncells
        phi = ones(1, ncells)*phi(1);
    end

    if length(theta) ~= ncells
        theta = ones(1,ncells)*theta(1);
    end

    for i = 1:ncells
        coords(:,:,i) = [com(i,:); com(i,:); com(i,:)];
        
        % phase offset
        offset = mod(offset,4);
        po = floor(offset(i))+1;
        
        % phase offset modulo
        pm = mod(offset(i)+1, po);
        
        if po > 1
            temp = W(1:(po-1),:,i);
            W(1:(5-po),:,i) = W(po:end,:,i);
            W((6-po):end,:,i) = temp;
            
            tempS = seq(1:(po-1),:,i);
            seq(1:(5-po),:,i) = seq(po:end,:,i);
            seq((6-po):end,:,i) = tempS;
        end

        % Determine arm lengths to keep mean arm length and allow cycle to 
        % start at the beginning of the defined sequence
        ts = first(seq(seq(:,1,i)~=0,1,i));        
        h1 = ((ts > 0)-0.5)*abs(ts)*deform;
        
        ls = first(seq(seq(:,2,i)~=0,2,i));
        g1 = ((ls > 0)-0.5)*abs(ls)*deform;
        h0(i) = L - h1;
        g0(i) = L - g1;
        if W(1,1,i) ~= 0
            h0(i) = h0(i) + pm*2*h1;
        else
            g0(i) = g0(i) + pm*2*g1;
        end
        
        % Apply to swimmer coordinates
        coords(1,1,i) = com(i,1)-h0(i); coords(3,1,i) = com(i,1)+g0(i);

        % Rotate to given angles
        coords(:,:,i) = racm3(coords(:,:,i), theta(i), phi(i));
        if fiber && ~heightspec
            coords(:,3,i) = coords(:,3,i) + r*1e4;
        end
    end

    if recenter
        coords = coords - mean(coords(2,:,:),3) + [0 0 mean(coords(2,3,:),3)];
    end
    
    % Arm bounds
    bounds = 0.5*deform*seq + L;

    %% Initialize
    % Number of steps
    nstep = round(ncycle*cyctime/dt);
    % nstep = ncycle*round(2*pi/omega / dt);

    % Arm lengths
    Lvec = zeros(nstep+1, 2, ncells);

    % Internal angle
    psivec = zeros(nstep+1, ncells);

    % Polar angle
    thetavec = zeros(nstep+1, ncells);

    % Azimuthal angle
    phivec = zeros(nstep+1, ncells);

    % Forces
    F = zeros(3*nstep,3, ncells);
    Fadh = zeros(3*nstep,3, ncells);

    % Velocities
    V = zeros(3*nstep,3,ncells);

    % Time & cycle
    tvec = zeros(nstep+1,1);
    cycvec = zeros(nstep+1,1);

    % All configurations
    allconf = zeros(3*(nstep+1),3,ncells);

    % First element
    psivec(1,:) = pi;
    thetavec(1,:) = theta; phivec(1,:) = phi;
    allconf(1:3,:,:) = coords;
    Lvec(1,:,:) = permute([h0 g0], [3 2 1]);
    econf = zeros(ncycle+1, 7, ncells);
    econf(1,2:4,:) = coords(2,:,:);
    econf(1,5:7,:) = permute([theta;phi;pi*ones(1,ncells)], [3 1 2]);


    % Phase counter
    phase = ones(ncells,1);

    % Initial xi distribution, W, deforming arm index & motion
    thisxi = zeros(ncells, 3);
    deforming = zeros(ncells, 1);
    motion = zeros(ncells, 1);
    for n = 1:ncells
        if all(abs(coords(:,3,n)-r) > 1e-3) && ~fiber
            thisxi(n,:) = 0;
        else
            thisxi(n,:) = changeXi(W(phase(n),:),xi,xibounds);
        end
        
        % Deforming arm
        deforming(n) = find(seq(1,:,n) ~= 0);
        
        % Motion
        motion(n) = 2*(any(seq(1,:,n) > 0) - 0.5);
    end
    
    % Initial W
    thisW = W(1,:,:);

    % Threshold hit
    hitthresh = false;
    
    
    arms = zeros(1,2,ncells);
    
    % Cycle counter
    k = zeros(ncells,1) + 2;
    
    % Phase switch boolean
    switchphase = false(ncells,1);
    

    i = 2;
    rdt = 0;
    t = 0;
    
    %% Simulation start
    while ~all(k > ncycle+1)        
        % Calculate velocities & forces
        [v, f, theta, phi, scale, fadh] = ...
            calc_vf(coords, eps, eta, thisW, thisxi, fthresh);
        
        % Threshold boolean
        if scale < 1
            hitthresh = true;
        end

        % Scale the timestep to speed up simulation
        useddt_temp = dt/scale;
        
        % Readout cap
        if rdt + useddt_temp > Rdt
            useddt_temp = (Rdt - rdt);
        end
        
        % Configuration update
        newcoords = coords + v*useddt_temp;

        % Arm lengths & Internal angle
        arms(1,1,:) = vecnorm(newcoords(2,:,:)-newcoords(1,:,:));
        arms(1,2,:) = vecnorm(newcoords(3,:,:)-newcoords(2,:,:));
            
        useddt = useddt_temp;
        
        for n = 1:ncells
            thisD = deforming(n);
            if motion(n)*(bounds(phase(n),thisD,n)-arms(1,thisD,n)) <= tol
                switchphase(n) = true;
                
                if motion(n)*(bounds(phase(n),thisD,n)-arms(1,thisD,n)) < -tol
                    % Smaller time step
                    dt_scale = (bounds(phase(n),thisD,n) - arms0(1,thisD,n)) /...
                        (arms(1,thisD,n) - arms0(1,thisD,n));
                    
                    useddt = useddt_temp*dt_scale;
                    % Redo with smaller step
                    newcoords = coords + v*useddt;
                    
                    if n > 1
                        switchphase(1:(n-1)) = false;
                    end
                    % Recalculate arms
                    arms(1,1,:) = vecnorm(newcoords(2,:,:)-newcoords(1,:,:));
                    arms(1,2,:) = vecnorm(newcoords(3,:,:)-newcoords(2,:,:));
                end
            end
        end
        
        arms0 = arms;
        coords = newcoords;
        
        t = t + useddt;
        rdt = rdt + useddt;
        
        % Save into readout matrices
        if rdt >= Rdt
            rdt = 0;
            F((3*i-2):(3*i),:,:) = f;
            V((3*i-2):(3*i),:,:) = v;
            Fadh((3*i-2):(3*i),:,:) = fadh;
            allconf((3*i-2):(3*i),:,:) = coords;
            tvec(i) = t;
            thetavec(i,:) = theta;
            phivec(i,:) = phi;
            psivec(i,:) = real(permute(acos(...
                dot(coords(1,:,:)-coords(2,:,:),...
                coords(3,:,:)-coords(2,:,:)) ./ ...
                (vecnorm(coords(1,:,:)-coords(2,:,:)).*...
                vecnorm(coords(3,:,:)-coords(2,:,:))...
                )), [2 3 1]));
            Lvec(i,:,:) = arms;
            
            i = i+1;
        end
        
        switched = find(switchphase);
        for j = 1:length(switched)
            if isempty(n)
                break
            end
            n = switched(j);
            phase(n) = mod(phase(n),4) + 1;
            
            % Update xi
            if all(abs(coords(:,3,n)-r) > 1e-3) && ~fiber
                thisxi(n,:) = 0;
            else
                thisxi(n,:) = changeXi(W(phase(n),:,n),xi,xibounds);
            end
            
            % Update W
            thisW(:,:,n) = W(phase(n),:,n);
            
            if phase(n) == 1
                econf(k(n),:,n) = [t ...
                    coords(2,:,n) ...
                    theta(n) phi(n) psivec(i-1,n)];
                k(n) = k(n)+1;
            end
            
            % Change deformer
            deforming(n) = find(seq(phase(n),:,n) ~= 0);
            % Motion
            motion(n) = 2*(any(seq(phase(n),:,n) > 0) - 0.5);
            % switchphase
            switchphase(n) = false;
            
        end

        % Collision w wall
        if any(coords(:,3) < collide)
            allconf = allconf(1:(i+1),:,:);
            tvec = tvec(1:(i+1));
            Lvec = Lvec(1:(i+1),:,:);
            psivec = psivec(1:(i+1),:);
            phivec = phivec(1:(i+1),:);
            thetavec = thetavec(1:(i+1),:);
            F = F(1:(3*i),:,:);
            Fadh = Fadh(1:(3*i),:,:);
            V = V(1:(3*i),:,:);
            econf = econf(1:(k-1),:,:);
            fprintf('Collision on step %1.0f\n', i)
            break
        end
    end
    econf = econf(1:(k-1),:,:);
    
    F = F(1:(3*(i-1)),:,:);
    V = V(1:(3*(i-1)),:,:);
    Fadh = Fadh(1:(3*(i-1)),:,:);
    allconf = allconf(1:(3*(i-1)),:,:);
    tvec = tvec(1:(i-1));
    thetavec = thetavec(1:(i-1),:);
    phivec = phivec(1:(i-1),:);
    psivec = psivec(1:(i-1),:);
    Lvec = Lvec(1:(i-1),:,:);
    % Subtract height if on fiber
    if fiber && ~heightspec
        econf(:,4,:) = econf(:,4,:) - r*1e4;
        allconf(:,3,:) = allconf(:,3,:) - r*1e4;
    end
    % Simulation end
end

% Extract first element (useful for expressions)
function x = first(X)
    x = X(1);
end