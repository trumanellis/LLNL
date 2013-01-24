Mixed Finite Element Lagrangian Hydrocode Testbed in Matlab

written by: 
Truman Ellis
ellis35@llnl.gov
ellis.truman@gmail.com

under supervision by: Robert Rieben (WCI) and Tzanio Kolev (CASC)

This work performed under the auspices of the U.S. Department of Energy by
Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344

INTRODUCTION
    This Lagrangian Hydrocode FEM testbed was written to explore the benefits 
of using higher order mixed finite elements in a lagrangian hydrocode.
    It has been written in a very modular way such that each problem script, 
such as Quad_Noh.m calls the same blocks of code for calculations thus 
reducing code repetition. I have tried to make the interface as simple as 
possible. Most of the important switches and definitions are defined at the 
top of each command script. 

SWITCHES AND DEFINITIONS
SaveFigures
    Pretty self-explanatory. If this is turned on, a subdirectory is created 
under FigureFiles and figures are saved according to the plotcycle, otherwise 
they just appear on the screen.

plotcycle
    Also pretty self-explanatory. This determines the frequency of figure 
plotting. Figure plotting takes time, so a higher number will be faster, but 
you will see less of the time history.

dtInit
    This is the initial timestep. It gets ramped up to a calculated stable
time step or dtMax, whichever is smaller.

dtMax
    This is the maximum allowable time step.

Method
    This defines the mixed FEM pair to use. The de facto standard in most 
current hydro codes is the Q1Q0 MFEM pair, in which velocity is defined in
the Q1 space at the nodes and the thermodynamic variables of pressure, 
density, and energy are defined in the Q0 space as constant values within 
each cell and discontinuous between cells.
    A higher order method that shows some promise is what we call the Q2Q1d
MFEM pair. The velocity is defined as a quadratice function with 9 dofs per cell.
The pressure and density are defined to be Q1 discontinuous which means that 
there are 4 dofs per cell, such that the pressure and density vary bilinearly 
within each cell, but are not discontinuous between cells. The energy remains
in the Q0 space. This finite element pair is similar to Taylor-Hood elements, 
but the discontinuous thermodynamic variables are more conducive to a large-
scale hydrocode with multi-materials.

isFullMassMatrixSolve
    isFullMassMatrixSolve decides whether to assemble a global mass matrix to
solve every cycle or to use mass lumping. The massQuadOrder must be set
appropriately. 
For mass lumping: Q1Q0 = 2, Q2Q1d = 3
For global mass matrix solve use any higher quadrature order

isMassUpdateEveryCycle
    isMassUpdateEveryCycle assumes mass lumping but reallocates nodal mass
based on the distortion of the elements

massQuadOrder
    The quadrature order to use for mass matrix or mass vector assembly

stiffQuadOrder
    The quadrature order to use for pressure gradient driven forces and 
artificial viscosity. A high enough stiffQuadOrder must be used to adequately 
capture the RHS.

Qfrac
    Use this to turn artificial viscosity on/off

hgfrac
    Use this to turn anti-hourglass forces on/off (Only works with Q1Q0)

NZx and NZy
    The number of zone in the x-direction and y-direction, respectively

xmin, xmax, ymin, ymax
    The mesh size

tstart and tstop
    the start and stop time for the simulation

jitter
    A random "jittering" factor -- this determines the magnitude of 
mesh distortion. This will usually be set to zero

maxcycle
    The maximum number of time steps to use

qquad and qlin
    These numbers are the coefficients for the artificial viscosity

DATA STRUCTURES
    Most data structures are defined such that the most basic elements of the
structure are the first index and more global features are later. For example, 
in cornerForce(i,j,k):
i refers to the x or y component of force
j refers to the location within the cell that the force is applied
k refers to the cell numbering in a global sense


