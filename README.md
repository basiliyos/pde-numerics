# pde_numerics
Spectral/pseudospectral PDE solvers, degenerate Schrödinger equation wavefront analysis

These MATLAB functions are intended for numerical simulation of PDEs in one space dimension, with emphasis placed on a degenerate nonlinear Schrödinger equation.  The PDEs are solved spectrally (or pseudospectrally in the case of nonlinear PDE).

The main script is DNLS.m.  This script calls other functions as needed, and each function call is in its own code section.  The first section of DNLS.m initializes many of the necessary variables and should always be run first.  Following that, the sections up to and including the function call PlayMovie.m are intended to be run sequentially.  After this, each function call requires multiple simulations of the PDE, and thus these functions may be called after only running the variable initialization section of DNLS.m.

# File functions and dependencies

DNLS.m: main script.  No input arguments or prerequisites.

  Setup.m: initializes spatial grid and initial data.  Can be run immediately after variable declarations.  If self-similar initial data is required, Setup.m calls function SelfSimInitialData.m.

    SelfSimInitialData.m: Solves ODE for self-similar solutions to DNLS, then uses solution to generate appropriate choice of initial data for original PDE.  Input parameter alpha determines choice of similarity variable.  Calls SimRHS.m when running ODE solver. 

      SimRHS.m: Right-hand side of self-similar ODE.
    
  SolvePDE.m: Solves choice of PDE for given initial data.  Can be run after Setup.m.  Acceptable choices of PDE can be found in PDEs_struct at beginning of DNLS.m.  Calls RHS.m when running ODE solver, and calls truncate.m if PDE is nonlinear.
  
    truncate.m: Truncates highest Fourier modes for pseudospectral method.  Requires degree of nonlinearity as input.
    
    RHS.m: Contains right-hand sides of various PDEs.  Calls deriv.m, and calls truncate.m for nonlinear PDE.
      
      deriv.m: Computes spectral derivative of a given order.
  
  Lpnorm.m: Computes and plots Lp norm of u at each time step.  Can be run after SolvePDE.m.

  WaterfallPlots.m: Plots waterfall plots of real part, imaginary part, and modulus of |u|^2 u_x (for DNLS solution).  Can be run after SolvePDE.m and after ux is computed (in section of DNLS.m immediately preceding Lpnorm.m call).

  Snapshots.m: Plots snapshots of modulus and phase of solution to PDE at various time intervals.  Can be run after SolvePDE.m.

  MoreSnapshots.m: Plots snapshots of other quantities relevant in DNLS analysis.  Calls deriv.m.  Can be run after Snapshots.m.

  WavefrontPosition.m: Computes and plots approximate location of wavefront as a function of time.  Can be calculated using either spatial or spectral method.  Includes option to omit plots.  

  WavefrontPlots.m: Plot real part, imaginary party, and modulus of |u|^2 u_x evaluated at wavefront as a function of time, as well as wavefront speed.  Can be run after WavefrontPosition.m.  Calls function timederiv.m to compute wavefront speed.

    timederiv.m: Computes time derivative of a function u on time grid tspan using Toeplitz differentiation matrix.
  
  SaveMovie.m: Saves movie of PDE solution to a given filename.  Can be run after SolvePDE.m.

  PlayMovie.m: Loops movie of PDE solution a prescribed number of times.  Can be run after SolvePDE.m.

  InitialDataStability.m: Examines stability of PDE with respect to variation in initial data.  Can be run after initial variable declarations in DNLS.m.  Calls functions Setup, SolvePDE, Lpnorm.

  LpResolutionStudy.m: Plots L^p norm comparison for PDE solution at different spatial resolutions.  Can be run after initiable variable declarations.  Calls Setup, SolvePDE, Lpnorm.
  
  AmpWidth: Plots wavefront position for various choices of amplitude and width of initial data.  Can be run after variable declarations.  Calls Setup, SolvePDE, WavefrontPosition.
  
  ViscConvergence: Plots wavefront position as a function of viscosity parameter.  Can be run after variable declarations.  Calls Setup, SolvePDE, WavefrontPosition.
  
  WavefrontResolutionStudy: Plots wavefront position at different spatial resolutions.  Can be run after variable declarations.  Calls Setup, SolvePDE, WavefrontPosition.
  
  WavefrontPowerlaw: Plots wavefront position for different choices of powerlaw initial data.  Can be run after variable declarations, but requires that u0switch = 'powerlaw'.  Calls functions Setup, SolvePDE, WavefrontPosition.

