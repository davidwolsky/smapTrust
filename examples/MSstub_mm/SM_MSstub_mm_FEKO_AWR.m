% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = 1e9;
fmax = 2e9;  
Nm = 41;	% Number of frequencies

% Initial input parameters 
xinit = [67.5]';  % Initial input parameters (only ls in this case)
xpinit = [2.1]';   % Initial implicit parameters (only eps_r in this case)

% filename = mfilename('SM_BandpassFilter.m')
filename = mfilename('.m');
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '')

% ------ Set up fine model ------ 
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'FEKO\']
%Mf.path = '/home/rib/Documents/masters/smap/smap/examples/MSstub/FEKO/';
Mf.name = 'MSstubOpen_mm';
Mf.solver = 'FEKO';
Mf.params = {'ls'};
% Mf.ximin = [0]';
% Mf.ximax = [1]';
Mf.ximin = [50]';
Mf.ximax = [90]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% ------  Set up coarse model (AWR) ------ 
Mc.path = [currentPath,'AWR\'];
Mc.name = 'MSstubCoarse_mm';
Mc.solver = 'AWR';
Mc.params = {'ls'};
Mc.Iparams = {'eps_r'};
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.xpmin = [1.8]';
Mc.xpmax = [3.0]';
Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 3;
% OPTopts.TRNi = OPTopts.Ni*2;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S1,1'};
OPTopts.globOpt = 0;
% OPTopts.globOpt = 1;
OPTopts.globOptSM = 0;
% OPTopts.globOptSM = 1;
OPTopts.globalSolver = 'ga';
OPTopts.optsGlobalOptim = optimoptions('ga');
OPTopts.optsGlobalOptim.Display = 'final';
OPTopts.localSolver = 'fmincon';
OPTopts.optsLocalOptim = optimoptions('fmincon');
OPTopts.optsLocalOptim.Display = 'iter-detailed';
% OPTopts.optsLocalOptim.DiffMinChange = 1e-4;
OPTopts.optsLocalOptim.DiffMinChange = 1e-5;
% OPTopts.optsLocalOptim.Diagnostics = 'on';
%
% OPTopts.goalType = {'minimax'};
% OPTopts.goalResType = {'S1,1_complex'};
% OPTopts.goalVal = {0.1-1j*0.2};
% OPTopts.goalWeight = {1};
% OPTopts.goalStart = {1.30e9};
% OPTopts.goalStop = {1.45e9};
% OPTopts.errNorm = {1};

% OPTopts.goalType = {'lt'};
% OPTopts.goalResType = {'S1,1_dB'};
% OPTopts.goalVal = {-20};
% OPTopts.goalWeight = {1};
% OPTopts.goalStart = {1.30e9};
% OPTopts.goalStop = {1.45e9};
% OPTopts.errNorm = {1};

% OPTopts.goalType = {'gt'};
% OPTopts.goalResType = {'S1,1_dB'};
% OPTopts.goalVal = {-10};
% OPTopts.goalWeight = {1};
% OPTopts.goalStart = {1.60e9};
% OPTopts.goalStop = {1.75e9};
% OPTopts.errNorm = {1};

OPTopts.goalType = {'lt', 'gt'};
OPTopts.goalResType = {'S1,1_dB', 'S1,1_dB'};
OPTopts.goalVal = {-20, -10};
OPTopts.goalWeight = {1, 0.1};
OPTopts.goalStart = {1.30e9, 1.60e9};
OPTopts.goalStop = {1.45e9, 1.75e9};
OPTopts.errNorm = {1,1};

OPTopts.plotIter = 1;
% Normalised tolerance in main loop
OPTopts.TolX = 0.1;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
% OPTopts.DeltaInit = 0.15;
OPTopts.DeltaInit = 0.25;
% OPTopts.DeltaInit = 0.35;
OPTopts.startWithIterationZero = 1;
% OPTopts.prepopulatedSpaceFile = '.mat';
OPTopts.verbosityLevel = 1;


% ------ Set up the SM ------ 
% The initial SM structure
Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA = 0;
SMopts.getB = 0;
SMopts.getc = 1;
SMopts.getG = 0;
SMopts.getxp = 1;
SMopts.getF = 0;
SMopts.getE = 0;
SMopts.getd = 0;

SMopts.globalSolver = 'ga';
SMopts.optsGlobalOptim = optimoptions('ga');
SMopts.optsGlobalOptim.Display = 'final';

SMopts.localSolver = 'fmincon';
SMopts.optsLocalOptim = optimoptions('fmincon');
SMopts.optsLocalOptim.Display = 'iter-detailed';
% SMopts.optsLocalOptim.Diagnostics = 'on';
% SMopts.normaliseAlignmentParameters = 0;
% SMopts.optsLocalOptim.DiffMinChange = 1e-5;
SMopts.normaliseAlignmentParameters = 1;
SMopts.optsLocalOptim.DiffMinChange = 1e-4;
% SMopts.optsLocalOptim.DiffMinChange = 2e-6;
% SMopts.plotAlignmentFlag = 0;
SMopts.plotAlignmentFlag = 1;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
% SMopts.errNorm = 1;
% An L2 norm is required here else the model steps in the wrong direction.
SMopts.errNorm = 2;
% SMopts.errNorm = 1;
% errW = 1
errW = zeros(size(Mf.freq));
% errW = ones(size(Mf.freq)).*0.2;
errW(Mf.freq > OPTopts.goalStart{1} & Mf.freq < OPTopts.goalStop{1}) = 1;
% errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 1
errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 0.1
SMopts.errW = errW;
% SMopts.wk = 5;
% SMopts.wk = 0; % <-
SMopts.wk = []

%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;

%% Notes:
% --- Alignment (PE) ---
% - If this model runs the wrong way it could be a tollerance issue with one of the optimiser options.
%   The minimum finite-difference gradient change may overstep initially. 
%   Alignment may not be working at all. (Shouldn't do this but could change the SMopts.errNorm...)
%   This is applicable to all of the paramets so if one is overstepped then it is a problem.
% - If there is no step away from the coarse model then the finite-difference gradient may be set too small.
%
% - It could also walk the wrong way if the goal weighting is not optimal. 
%   This is especially true for goals that initially sit in a null.
