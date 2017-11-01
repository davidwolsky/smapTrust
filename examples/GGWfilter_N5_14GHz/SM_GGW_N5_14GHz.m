% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = 13.8e9;
fmax = 14.2e9;  
Nm = 201;	% Number of frequencies

% Initial input parameters
xinit = [0, 0, 0, 1, 1, 1]';  
% Initial implicit parameters
xpinit = [0,1,1,1,1,1,1,-1.04,-1.04,-1.04]';
% xpinit = [-0.5953    0.7130    0.9302    1.0320    1.1093    0.8660    0.9690   -0.8992   -1.0866 -1.0925]';

filename = mfilename('SM_GGW_N3_14GHz.m');
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '');

% Set up fine model
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'CST\'];
Mf.name = 'GGW_N5_14GHz';
Mf.solver = 'CST';
Mf.params = {'L1_term', 'L2_term', 'L3_term', 's12_fact', 's23_fact', 'x_fact'};
Mf.ximin =  [-0.5, -0.25, -0.25, 0.3, 0.5, 0.5]';
Mf.ximax =  [0.5,   0.25,  0.25, 1.5, 1.5, 2]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (AWR)
Mc.path = [currentPath,'AWR\'];
Mc.name = 'GGW_N5_14GHz';
Mc.solver = 'AWR';
Mc.params = Mf.params;
Mc.ximin =  Mf.ximin;
Mc.ximax =  Mf.ximax;
Mc.Iparams ={'P','as_fact','bs_fact','cs_fact','ax_fact','bx_fact','cx_fact','L1c','L2c','L3c'};
% Mc.xpmin =  [-0.9,0.7,0.8,0.9,0.85,0.8,0.85,-1.1,-1.15,-1.15]';
% Mc.xpmax =  [ 0.9,1.1,1.1,1.2,1.15,1.1,1.15,-0.8, -1,  -1]';
Mc.xpmin =  [-0.9,0.7,0.8,0.9,0.85,0.8,0.85,-1.1,-1.2,-1.2]';
Mc.xpmax =  [ 0.9,1.1,1.1,1.2,1.15,1.1,1.15,-0.3,-1,-1]';
Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up the SM
% The initial SM structure
% Sinit.c = [-0.72,-1.22,0,0];
% Sinit = [];
Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA = 0;
SMopts.getB = 0;
SMopts.getc = 0;
SMopts.getG = 0;
SMopts.getxp = 1;
SMopts.getF = 0;
SMopts.getE = 0;
SMopts.getd = 0;

%% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 5;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S1,1'};
OPTopts.globOpt = 1;

OPTopts.globalSolver = 'ga';
OPTopts.optsGlobalOptim = optimoptions('ga');
OPTopts.optsGlobalOptim.Display = 'final';
OPTopts.optsGlobalOptim.Display = 'iter';
OPTopts.optsGlobalOptim.PopulationSize = 50;
OPTopts.optsGlobalOptim.Generations = 6;

% OPTopts.localSolver = 'fmincon';
% OPTopts.optsLocalOptim = optimoptions('fmincon');
% OPTopts.optsLocalOptim.Display = 'iter-detailed';
% OPTopts.optsLocalOptim.Diagnostics = 'on';
% OPTopts.optsLocalOptim.DiffMinChange = 1e-5;

OPTopts.localSolver = 'fminsearchcon';
OPTopts.optsLocalOptim.Display = 'iter';


OPTopts.goalType =      {'gt', 'lt', 'gt'};
OPTopts.goalResType =   {'S1,1_dB','S1,1_dB','S1,1_dB'};
OPTopts.goalVal =       {-10, -16.4, -10};
OPTopts.goalWeight =    {1,1,1};
OPTopts.goalStart =     {fmin, 13.93e9, 14.08e9};
OPTopts.goalStop =      {13.92e9, 14.07e9, fmax};
OPTopts.errNorm =       {1, 1, 1};

OPTopts.plotIter = 1;
OPTopts.TolX = 1e-2;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
% OPTopts.alp2 = 0.75;
OPTopts.DeltaInit = 0.25;
OPTopts.startWithIterationZero = true;
% OPTopts.prepopulatedSpaceFile = 'SMLog_bandpassFilter.mat'

%% Set up the parameter extraction
OPTopts.globOptSM = 1;

% SM global optimiser
SMopts.globalSolver = 'ga';
SMopts.optsGlobalOptim = optimoptions('ga');
SMopts.optsGlobalOptim.Display = 'final';
% gaOptions = gaoptimset('PopulationSize',100,'Generations',2);
SMopts.optsGlobalOptim.Display = 'iter';
SMopts.optsGlobalOptim.PopulationSize = 200;
SMopts.optsGlobalOptim.Generations = 10;

% SM local optimiser

% SMopts.localSolver = 'fmincon';
% SMopts.optsLocalOptim = optimoptions('fmincon');
% % SMopts.optsLocalOptim = optimoptions('fmincon','Algorithm','sqp');
% SMopts.optsLocalOptim.Display = 'iter-detailed';
% SMopts.optsLocalOptim.Diagnostics = 'on';
% SMopts.normaliseAlignmentParameters = 1;
% SMopts.optsLocalOptim.DiffMinChange = 1e-5;

SMopts.localSolver = 'fminsearchcon';
SMopts.optsLocalOptim.Display = 'iter';

% SMopts.localSolver = 'fminsearch';
% SMopts.optsLocalOptim.Display = 'iter';

SMopts.plotAlignmentFlag = 1;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.errNorm = 1;
% errW = 1;
errW = ones(size(Mf.freq)).*0.3;
% errW(Mf.freq > OPTopts.goalStart{1} & Mf.freq < OPTopts.goalStop{1}) = 1
errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 1;
% errW(Mf.freq > OPTopts.goalStart{3} & Mf.freq < OPTopts.goalStop{3}) = 1
SMopts.errW = errW;
% SMopts.wk = 5;
SMopts.wk = []; % <---
% SMopts.wk = 1;


%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;

