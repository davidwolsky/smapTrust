% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = 5e9;
fmax = 20e9;  
Nm = 151;	% Number of frequencies

% Initial input parameters 
%           l1      l2      s
xinit = [   66.727, 60.228, 9.592]';
% xinit = [   78.964, 81.210, 7.901]';
%           cm
xpinit = [  0.0445, 30]'; 

% filename = mfilename('SM_BandpassFilter.m')
filename = mfilename('.m');
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '')

% ------ Set up fine model ------ 
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'FEKO\']
%Mf.path = '/home/rib/Documents/masters/smap/smap/examples/MSstub/FEKO/';
Mf.name = 'DoubleFoldedStub';
Mf.solver = 'FEKO';
Mf.params = {'l1',  'l2',   's'};
Mf.ximin = [ 25.0,   25.0    1.0]';
Mf.ximax = [ 90.0    90.0    15.0]';

%%%%%%%%%%%% TODO_DWW: Try this next

% Mf.ximax = [ 90.0    90.0    30.0]';
Mf.freq = reshape(linspace(fmin,fmax,Nm), Nm, 1);

% ------  Set up coarse model (AWR) ------ 
Mc.path = [currentPath,'AWR\'];
Mc.name = 'DoubleFoldedStub';
Mc.solver = 'AWR';
Mc.params = {'l1',  'l2',   's'};
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.Iparams = {'cm',     'l3'};
Mc.xpmin =   [0.025,     20.0]';
Mc.xpmax =   [0.080,     40.0]';
Mc.freq = reshape(linspace(fmin,fmax,Nm), Nm, 1);

% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 3;
% OPTopts.TRNi = OPTopts.Ni*2;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S2,1'};
OPTopts.globOpt = 0;
% OPTopts.globOpt = 1;
OPTopts.globOptSM = 0;
% OPTopts.globOptSM = 1;
OPTopts.globalSolver = 'ga';
OPTopts.optsGlobalOptim = optimoptions('ga');
OPTopts.optsGlobalOptim.Display = 'final';
OPTopts.localSolver = 'fminsearchcon';
% OPTopts.optsLocalOptim = optimoptions('fmincon');
OPTopts.optsLocalOptim.Display = 'iter';
% OPTopts.optsLocalOptim.DiffMinChange = 1e-4;
% OPTopts.optsLocalOptim.DiffMinChange = 1e-5;
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

OPTopts.goalType =      {'gt',      'lt',       'gt'};
OPTopts.goalResType =   {'S2,1_dB', 'S2,1_dB',  'S2,1_dB'};
OPTopts.goalVal =       {-3,        -30,        -3};
OPTopts.goalWeight =    {1,         1,          0.1};
OPTopts.goalStart =     {5.0e9,     12.0e9,     16.5e9};
OPTopts.goalStop =      {9.5e9,     14.0e9,     20.0e9};
OPTopts.errNorm =       {1,         1,          1};

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
% OPTopts.prepopulatedSpaceFile = 'SMLogMSstubOpen_mm_EXAMPLE.mat';
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

SMopts.localSolver = 'fminsearchcon';
% SMopts.optsLocalOptim = optimoptions('fminsearchcon');
SMopts.optsLocalOptim.Display = 'iter';
% SMopts.optsLocalOptim.Diagnostics = 'on';
% SMopts.normaliseAlignmentParameters = 0;
% SMopts.optsLocalOptim.DiffMinChange = 1e-5;
SMopts.normaliseAlignmentParameters = 1;
% SMopts.optsLocalOptim.DiffMinChange = 1e-4;
% SMopts.optsLocalOptim.DiffMinChange = 2e-6;
% SMopts.plotAlignmentFlag = 0;
% SMopts.plotAlignmentFlag = 1;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
% SMopts.errNorm = 1;
% An L2 norm is required here else the model steps in the wrong direction.
% SMopts.errNorm = 2;
% SMopts.errNorm = 1;
errW = zeros(size(Mf.freq));        % --- use this ---
errW(Mf.freq > OPTopts.goalStop{1} & Mf.freq < OPTopts.goalStart{3}) = 10;
% errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{3}) = 1;
errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 20;
errW(Mf.freq > 12.5e9 & Mf.freq < 13.5e9) = 50;
SMopts.errW = errW;                 % --- use this ---
% SMopts.wk = 5;
SMopts.wk = 0; % <-
% SMopts.wk = 3;
% SMopts.wk = []


plotOpts = {};
plotOpts.plotModelOpts = {};
plotOpts.plotModelOpts.ylim = [-80 0];
plotOpts.plotModelOpts.pbaspectVec = [1 1 1];


%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit, Sinit, SMopts, Mf, Mc, OPTopts, plotOpts);

keyboard;