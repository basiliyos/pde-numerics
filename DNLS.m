% Main script: call other functions to set up and simulate PDEs, perform
% wavefront analysis, L^p norm calculations, resolution and convergence
% studies, etc.

%% Define all relevant variables and structures

% Keys for PDE type struct
keys = {'exponential'; 'characteristic'; 'heat equation'; ... 
    'free linear Schrodinger'; 'linear Schrodinger with potential'; ...
    'KdV'; 'Burgers'; 'cubic focusing NLS'; 'cubic defocusing NLS'; ... 
    'modified DNLS'; 'DNLS'; 'DNLS pad'; 'DNLS parabolic regularization'};

% Values for PDE type struct (for degree of nonlinearity of each PDE)
values = [1; 1; 1; 1; 1; 2; 2; 3; 3; 3; 3; 3; 3];

% Form PDE struct
PDEs_struct = containers.Map(keys, values);

% Integer to keep track of current figure number
i = 1;


% Basis for spatial resolution
k = 9;


% Timestep
dt = 0.0025;


% Maximum time
tmax = 0.3;


% Type of initial data
u0switch = 'parabola';


% Initial data amplitude
amp0 = 1;


% Initial data width
width = pi/2;


% Power for powerlaw initial data
pow = 1/3;


% Type of ODE solver to use
solver = '23tb';


% Absolute tolerance for ODE solver options
abstol = 0;


% Relative tolerance for ODE solver options
reltol = 0;


% Name of PDE to solve
PDE_name = 'DNLS';


% Amount of numerical viscosity to use (on the order of spatial stepsize *
% initial amplitude)
visc = 4;


% Get nonlinearity degree from PDE struct
deg_nonlinearity = PDEs_struct(PDE_name);


% Function calls requiring multiple initial functions:

% First function: type, amplitude, width, power
firstu0switch = 'parabola';
firstamp0 = 1;
firstwidth = pi/2;
firstpow = 0;

% Second function: type, amplitude, width, power
secondu0switch = 'sim';
secondamp0 = 1;
secondwidth = pi/2;
secondpow = 0;


% Choice of p in L^p norms
p = 2;

% Number of snapshots to plot 
num_snapshots = 5;

% Cutoff for computing phi near wavefront
phicutoff = 1/20;

% Relative spatial cutoff for wavefront location determination
relativespatialcutoff = 1000;


% Relative spectral cutoff for wavefront location determination
relativespectralcutoff = 3;


% Number of times to loop movie
loops = 10;


% Filename to save movie to
filename = 'DNLS_Solution.avi';


% Resolution vector for resolution study
kvec = [8:10];


% Amplitude vector for wavefront study
ampvec = [1,3/2];


% Width vector for wavefront study
widthvec = [1/2];


% Choice of wavefront calculation method for amp/width wavefront study and
% viscosity wavefront study
method = 'spectral';

% Choose appropriate cutoff based on method
if strcmp(method,'spectral')
    relativecutoff = relativespectralcutoff;
end
if strcmp(method,'spatial')
    relativecutoff = relativespatialcutoff;
end


% Vector of viscosity values (to study effect of viscosity on wavefront)
viscvec = [1,3,5,10];


% Fraction of indices to consider when performing linear fit 
minindexfrac = 0.2; 
maxindexfrac = 0.8;


% Vector of choices of power initial data
powvec = [1/6, 1/4, 1/3, 1/2];

%% Call setup function for grids, initial data
[N, h, xmin, xmax, x, tspan, u0, nplots, maximum] = ... 
    Setup(k, dt, tmax, u0switch, amp0, width, pow);



%% Call function to solve PDE
[t, u] = SolvePDE( u0, PDE_name, solver, abstol, reltol, ...
    x, tspan, deg_nonlinearity, N, h, amp0, visc );



%% Compute first spatial derivative of u
ux = zeros(size(u));
for j = 1:length(tspan)
    ux(j,:) = deriv(u(j,:)',N,1)';    
end



%% Compute and plot Lp norm of solution as a function of time
[u_p, i] = Lpnorm(u, x, tspan, p, i);



%% Waterfall plots
i = WaterfallPlots( u, ux, x, t, i );



%% Snapshot plots of modulus and phase
[r, phi, i] =  ... 
    Snapshots( u, x, xmin, xmax, maximum, t, num_snapshots, phicutoff, i );



%% More snapshot plots
[ rx, phix, i ] = MoreSnapshots( r, phi, x, t, num_snapshots, N, i );




%% Spatial cutoff method for finding wavefront location
[wavefrontspatial, wavefrontindexspatial, i] = ... 
    WavefrontPosition(u, N, t, x, xmin, xmax, relativespatialcutoff, 'spatial', phix, i, ... 
    amp0, width, 'yes');




%% Wavenumber method for finding wavefront location
[wavefrontspectral, wavefrontindexspectral, i] = ... 
    WavefrontPosition(u, N, t, x, xmin, xmax, relativespectralcutoff, 'spectral', phix, i, ... 
    amp0, width, 'yes');




%% Plot various quantities of u at wavefront
[ i ] =  WavefrontPlots( wavefrontspectral, wavefrontindexspectral, ... 
    u, ux, t, i );




%% Save movie to file
i = SaveMovie(u, x, xmin, xmax, maximum, t, filename, nplots, i);




%% Play movie showing Re(u), Im(u), and |u|
i = PlayMovie( u, x, xmin, xmax, maximum, t, nplots, loops, i );




%% Call function to compare PDE solutions for different initial data
[error, i] = InitialDataStability(k, dt, tmax, visc, ...
    PDE_name, solver, deg_nonlinearity, abstol, reltol, ...
    firstu0switch, firstamp0, firstwidth, firstpow, ... 
    secondu0switch, secondamp0, secondwidth, secondpow, i);




%% L^p norm resolution study and powerlaw fit
[i, power] = LpResolutionStudy(kvec, dt, tmax, u0switch, amp0, width, pow, ... 
    PDE_name, solver, abstol, reltol, deg_nonlinearity, visc, i);




%% Amplitude/width wavefront studies
[ i, amp_width_wavefront ] = AmpWidth( k, xmin, xmax, tmax, dt, pow, visc, ... 
    PDE_name, deg_nonlinearity, solver, abstol, reltol, u0switch, ... 
    relativecutoff, method, ampvec, widthvec, i);




%% Convergence of wavefront approximation as function of viscosity
[ i, visc_wavefront ] = ... 
    ViscConvergence( viscvec, k, dt, tmax, u0switch, relativecutoff, method, ... 
    amp0, width, pow, PDE_name, solver, abstol, reltol, deg_nonlinearity, i );




%% Wavefront resolution study
[ i, wavefront_resolution ] =  ... 
    WavefrontResolutionStudy(kvec, dt, tmax, u0switch, amp0, width, pow, ... 
    PDE_name, solver, abstol, reltol, deg_nonlinearity, visc, ... 
    relativecutoff, method, i );




%% Wavefront powerlaw study and linear fit analysis (return relativeerrorvec)
[ i, stufftoplot, relativeerrorvec ] = ... 
    WavefrontPowerlaw( powvec, dt, tmax, u0switch, amp0, width, k, ... 
    PDE_name, solver, abstol, reltol, deg_nonlinearity, visc, ... 
    relativecutoff, method, i, minindexfrac, maxindexfrac  );
