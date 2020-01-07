%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sphere3
Melissa Mai

Rigid three-sphere swimmer/crawler (no internal bending), with 
zero-z-velocity constraint when xi is active. Three-sphere swimmer based 
off Najafi & Golestanian (2004). Surface constraints based off 
Daddi-Moussa-Ider et al (2018). Regularized Blake tensors derived in
Ainley et al (2018).

This model uses discrete, constant-velocity deformation phases. Flow
field calculations are not included in this script.


INPUTS
'L', 'arm', 'length'      Mean length of each arm. Extensions and
                              contractions, in addition to initial
                              configuration, will be determined to
                              maintain a mean arm length of L.
                              Default: 1
'dt'                      Initial time step
                              Default: 1e-2
'r', 'a'                  Radius of the beads. Assumes the same radius
                              for all beads.
                              Default: 0.1*L
'eta', 'viscosity'        Viscosity
                              Default: 1
'h', 'height', 'z0',      Initial height of the swimmer's center of mass
  'heights', 'wall'           off the wall (which is assumed to be the 
                              xy-plane)
                              Default: r
'com', 'cm', 'init_com',  Initial center of mass. Can be given as a
  'initcom', 'config',        2-vector to specify projection onto the
  'initconfig',               xy-plane and use the 'height' parameter for
  'init_config'               the z-component, or can be given as a
                              3-vector to override the 'height'
                              parameter.
                              Default: [0 0 height]
'ncycle', 'ncycles'       Total number of swimming cycles
                              Default: 10
'xi', 'adh', 'adhesion'   Global adhesive friction parameter. Will be
  'friction'                  scaled at each phase according to the
                              movement:
                              Leading+    [high high low]
                              Leading-    [low high high]
                              Trailing+   [high high low]
                              Trailing-   [low high high]
                              Change 'xibounds', 'xihigh', and/or 'xilow'
                              to alter these scales
                              Default: 0
'u0', 'dv', 'w'           Deformation/distortion velocity for the arms.
                              Will use the same value of u0 for each
                              deformation; to change this, see 'seq'.
                              Default: 0.5*L
'deform'                  Length of full contraction/extension. For
                              example, a 'deform' of 0.5 with a mean arm
                              length L of 10 leads to arm lengths ranging
                              from 7.5 to 12.5.
                              Default: 0.5*L
'seq', 'sequence'         Sequence of arm motions. Typically given in
                              terms of +/- 1 but can be changed to
                              modulate the deformation velocity of each
                              step. Each row represents a phase, and the
                              columns represent [trailing leading] arms.
                              Use +1 for extension and -1 for
                              contraction. The initial configuration will
                              be determined so that the first step can be
                              completed.
                              Default: [1 0; # Trailing extension
                                        0 1; # Leading extension
                                       -1 0; # Trailing contraction
                                        0 -1]# Leading contraction
'phase', 'startphase',    Phase offset, in the range [0 4] (though if
  'phaseoffset',              given as a value greater than four, it will
  'offset'                    be rescaled via modulo to be within that
                              range. Can be given as a vector to give
                              individual offset values to each cell or as
                              a single value to assign a global phase
                              offset. The offset will determine the
                              starting configuration of the cell, with
                              [0 1) corresponding to the first phase,
                              [1 2) corresponding to the second phase,
                              etc. For example, an offset of 0.5 will
                              start the swimmer halfway through the first
                              phase defined by seq (see above). An offset
                              of 2 will begin the cell at the start of
                              the third motion phase.
                              Default: 0
'fthresh', 'thresh'       Threshold force. All forces will be scaled if
                              any of the required forces (the norms, not
                              the components) exceeds the threshold.
                              Forces are scaled as 
                                  F * min(1, fthresh/max(fnorms))
                              Default: 100
'theta', 'polar'          Polar angle (angle wrt z-axis). 
                              Default: pi/2
'phi', 'azimuthal',       Azimuthal angle (angle in xy-plane wrt to
  'azimuth'                   x-axis)
                              Default: 0
'collide', 'collision'    Collision threshold. If any sphere approaches
                              within 'collide' of the wall, the
                              simulation stops. A notification of the
                              collision and the step at which it occurs
                              will be generated in the command line.
                              Default: 0.75*r
'xihigh', 'high', 'xih'   Upper bound for the adhesion scaling (see
                              'xi'). Overridden if 'xibounds' is defined.
                              Default: 1
'xilow', 'low', 'xil'     Lower bound for the adhesion scaling (see
                              'xi'). Overridden if 'xibounds' is defined.
                              Default: 0.2
'xibounds'                Upper and lower bounds for adhesion scaling
                              (see 'xi'). Of the form [high low], though
                              it will be rearranged if [low high] is
                              given. Overrides 'xihigh' and 'xilow' if
                              given, otherwise, will be defined as
                              [xihigh xilow].
                              Default: [xihigh xilow]
'fiber'                   Boolean expression whether the simulation is
                              run on a fiber (ie, substrate hydrodynamics
                              are turned off by setting cell height
                              arbitrarily large (r*1e4), but adhesion
                              still exists. Height will be reset back to 
                              height at the end of the simulation.) If 
                              'height' parameter is explicitly defined, 
                              the cell will still be detached from the 
                              wall, but at the specified height (and not 
                              arbitrarily far away)
                              Default: false
'recenter'                Recenter the starting configuration so it is
                              centered at the xy origin
                              Default: false;
'tol', 'tolerance'        Tolerance for boundary cutoff; the motion phase
                              will advance when the arm length is within
                              tol of the actual value
                              Default: 1e-6
'rdt', 'readout'          Readout time step. The program will not record
                              every step (for outputs like allconf, V, F,
                              thetavec, phivec, psivec, Lvec, Fadh, tvec)
                              but will report in intervals of (at least)
                              rdt.
                              Default: 100*dt

NOTE: Input names are not case-sensitive.
NOTE: Presets of inputs can be defined and passed as a single cell: see
      example. Extra inputs can be added after the preset, but only one
      preset can be defined.


OUTPUTS
econf                     (ncycle+1)x7xncell matrix of the time, 
                              coordinates of the middle sphere, and 
                              angles at the end of each cycle. If the 
                              cell collides with the wall (see 'collide' 
                              above), this matrix will be truncated to 
                              only include completed cycles. The matrix 
                              includes the initial position at t=0.
                              [time x2 y2 z2 theta phi psi]
allconf                   (3*nstep+3)x3xncell matrix of the coordinates 
                              at each time step. Coordinates at each time
                              step are each triplet of rows (each row is
                              one sphere).
V                         (3*nstep+3)x3xncell matrix of the velocity 
                              components at each time step. Same
                              dimensions as allconf.
F                         (3*nstep+3)x3xncell matrix of the force 
                              components at each time step. Same
                              dimensions as allconf.
thetavec                  (nstep+1)xncell vector of the polar angles at 
                              each time step
phivec                    (nstep+1)xncell vector of azimuthal angles
psivec                    (nstep+1)xncell vector of internal angles
Lvec                      (nstep+1)x2xncell matrix of arm lengths. 
                              Columns given as [trailing leading]
tvec                      (nstep+1)x1 vector of time points
cycvec                    (nstep+1)x1 vector of time points transformed
                              to their relative position within the cycle
hitthresh                 Boolean of if the threshold force (fthresh) was
                              hit


USAGE
Ex: >> sphere3
  Will run the simulation with all default parameters

Ex: >> sphere3('L', 20, 'height', 4, 'phi', pi/2)
  Will use a mean arm length of 20 and will start the swimmer 4 units
  above the wall with azimuthal angle pi/2 (swimming in the +y direction)

Ex: >> preset1 = {'initconfig', [0 0 0; 1 3 0], 'phi', [0 pi], ... 
                  'ncycle', 200, 'dt', 1e-3, 'rdt', 1e-1}
    >> sphere3(preset1)
  Allows pre-defined sets of inputs to be passed as a single input. This
  example initializes two anti-parallel cells, one positioned at [0 0 0]
  moving in the +x direction and the other at [1 3 0] moving in the -x
  direction. Runs for 200 cycles with dt = 1e-3 and readout time step of
  1e-1. Additional inputs can be included after the preset; ie, 
      >> sphere3(preset1, 'L', 10)
  is permissible, but only one preset is allowed. 
      >> sphere3(preset1, preset2)
  is not allowed.


SEE ALSO
  - calc_vf
  - racm3
  - changeXi
  - build_vmat
  - Smat
  - vcurve
  - phicurve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


calc_vf

Calculates the forces and velocities for individual beads in the
three-sphere model using a set of defined constraints:

For swimming:
  Force-free (all three directions)
  Torque-free (beta- and gamma- components)
  Deformation velocities
  Rigid body

For crawling:
  Force-free (x & y)
  Torque-free (beta)
  Zero z-velocities
  Deformation velocities
  Rigid body

Then calculates velocities using V = M*F, with M the mobility matrix
(calcuated with build_vmat).


INPUTS
coords                    Coordinates for the cells, 3x3xncell
                              [x1 y1 z1;
                               x2 y2 z2;
                               x3 y3 z3]
eps                       Epsilon (regularization scale) / bead radius
eta                       Viscosity
thisW                     1x2xncell vector of the arm velocities
                              [trailing leading]
thisXi                    ncellx3 matrix of xi values
fthresh                   Threshold force


OUTPUTS
V                         3x3xncell matrix of velocity components
F                         3x3xncell matrix of internal force components
theta                     1xncell vector of polar angles
phi                       1xncell vector of azimuthal angles
scale                     Threshold scaling factor, 
                              min(1, fthresh / max(|F|))
Fadh                      3x3xncell matrix of adhesive force components


SEE ALSO
  - build_vmat
  - Smat
  - sphere3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


build_vmat

Build the velocity matrix for a three-sphere swimmer under the influence
of a planar surface's hydrodynamics. 
Hydrodynamics adapted from Ainley et al (2008)
Wall is assumed to be the xy-plane (z = 0)


INPUTS
coords                    Coordinates of all three spheres, 3x3xncell
                              [x1 y1 z1;
                               x2 y2 z2;
                               x3 y3 z3]
ncells                    Number of cells
eta                       Viscosity
eps                       Epsilon, blob cutoff parameter (r)
Xi                        1x(9*ncell) vector of xi values 
                              [xi1x xi1y xi1z ... xi3y xi3z]

OUTPUT
vmat                      Full (9*ncell)x(9*ncell) matrix,
                              [S1->1   S2->1  S3->1
                               S1->2   S2->2  S3->2
                               S1->3   S2->3  S3->3]

                          = [(s1x s1y s1z s2x s2y s2z s3x s3y s3z)->1x
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->1y
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->1z
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->2x
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->2y
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->2z
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->3x
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->3y
                             (s1x s1y s1z s2x s2y s2z s3x s3y s3z)->3z],

                              where Sm->n is the mobility submatrix 
                              describing the effect of particle m on n, 
                              and (smi)->nj is the individual component
                              describing the effect of particle m's 
                              i-directional force on particle n's 
                              j-directional velocity. Sn->n is the self 
                              mobility tensor.
vcube                     Deconstructed mobility matrix,
                              3x3x(3*ncell)x(3*ncell). vcube(:,:,n,m) is
                              S(m->n).


SEE ALSO
  - Smat
  - calc_vf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Smat

Calculates the mobility matrix for regularized Stokeslets with a 
wall boundary condition using the method of images, as described in 
Ainley et al (2008). Only for one pair. Use with build_vmat for full
system.
Wall is assumed to be the xy-plane (z = 0).


INPUTS
Xe                        Coordinates of point at which velocity will be 
                              evaluated
Xs                        Coordinates of the Stokeslet (force)
eps                       Epsilon, blob cutoff parameter
eta                       Viscosity


OUTPUT
S                         Mobility matrix for s -> e
                              ie, v(s->e) = S*fs

SEE ALSO
  - build_vmat
  - vfield_S


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


racm3

Rotate a 3-sphere swimmer/crawler about its center of mass, 3d

INPUTS
coords                    3x3 (row corresponding to point)
theta                     polar angle (wrt z-axis)
phi                       azimuthal angle (wrt to x-axis in xy-plane)


OUTPUT
newcoords                 3x3 matrix of coordinates for rotated cell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


changeXi

Updates the adhesion on each bead depending on the motion:
  Trailing arm extension:     [high high low]
  Leading arm extension:      [high high low]
  Trailing arm contraction:   [low high high]
  Leading arm contraction:    [low high high]

INPUTS
thisW                     Current motion [trailing leading]
xi                        Global adhesion parameter
xibounds                  [High Low] scale factors for adhesion


OUTPUT
thisXi                    1x3 vector of the adhesion values for the cell
