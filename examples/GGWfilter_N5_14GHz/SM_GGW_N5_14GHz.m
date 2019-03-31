% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = 13.8e9;
fmax = 14.2e9;  
Nm = 401;	% Number of frequencies

% Initial input parameters
xinit = [0, 0, 0, 0, 0, 0]';  
% Initial implicit parameters
xpinit = [0.5282]';

filename = mfilename('SM_GGW_N5_14GHz.m');
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '');

% Set up fine model
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'CST\'];
Mf.name = 'GGW_N5_14GHz';
Mf.solver = 'CST';
Mf.params = {'L1_term', 'L2_term', 'L3_term', 's12_term', 's23_term', 'x_term'};
Mf.ximin =  [-0.5, -0.5, -0.5, -0.5, -0.5, -1.0]';
Mf.ximax =  [ 0.5,  0.5,  0.5,  1.0,  1.0,  1.0]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (AWR)
% Mc.path = [currentPath,'AWR\'];
% Mc.name = 'GGW_N5_14GHz';
% Mc.solver = 'AWR';
% Mc.params = Mf.params;
Mc.path = [currentPath,'MATLAB\'];
addpath([Mc.path]);
Mc.name = @GGW_N5_14GHz;
Mc.solver = 'MATLAB';
Mc.params = {'x';'xp';'f'};
Mc.ximin =  Mf.ximin;
Mc.ximax =  Mf.ximax;
Mc.Iparams ={'P'};
Mc.xpmin =  [0]';
Mc.xpmax =  [1.8]';
Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up the SM
% The initial SM structure
Sinit.c = [0.0551   -0.0183   -0.0214    0.0417    0.0289   -0.2441]';
% Sinit = [];
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

%% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 15;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S1,1'};
OPTopts.globOpt = 0;
OPTopts.useScAsOpt = true;

OPTopts.globalSolver = 'ga';
OPTopts.optsGlobalOptim = optimoptions('ga');
OPTopts.optsGlobalOptim.Display = 'final';
OPTopts.optsGlobalOptim.Display = 'iter';
OPTopts.optsGlobalOptim.PopulationSize = 300;
OPTopts.optsGlobalOptim.Generations = 8;

% OPTopts.localSolver = 'fmincon';
% OPTopts.optsLocalOptim = optimoptions('fmincon');
% OPTopts.optsLocalOptim.Display = 'iter-detailed';
% OPTopts.optsLocalOptim.Diagnostics = 'on';
% OPTopts.optsLocalOptim.DiffMinChange = 1e-5;

OPTopts.localSolver = 'fminsearchcon';
OPTopts.optsLocalOptim.Display = 'iter';


OPTopts.goalType =      {'gt', 'lt', 'gt'};
OPTopts.goalResType =   {'S1,1_dB','S1,1_dB','S1,1_dB'};
OPTopts.goalVal =       {-5, -16.4, -5};
OPTopts.goalWeight =    {1,1,1};
OPTopts.goalStart =     {fmin, 13.93e9, 14.08e9};
OPTopts.goalStop =      {13.92e9, 14.07e9, fmax};
OPTopts.errNorm =       {1, 1, 1};

OPTopts.plotIter = 1;
OPTopts.TolX = 1e-3;

% Switch of the TR...
OPTopts.TRNi = 1;
OPTopts.eta1 = 0.05*0;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 1;%2.5;
OPTopts.alp2 = 1;%0.25;
OPTopts.DeltaInit = 1;%0.25;
OPTopts.startWithIterationZero = true;
% OPTopts.prepopulatedSpaceFile = 'SMLogGGW_N5_14GHz_20171218_2328.mat';

%% Set up the parameter extraction
OPTopts.globOptSM = 2;

% SM global optimiser
% SMopts.globalSolver = 'ga';
% SMopts.optsGlobalOptim = optimoptions('ga');
% SMopts.optsGlobalOptim.Display = 'final';
% % gaOptions = gaoptimset('PopulationSize',100,'Generations',2);
% SMopts.optsGlobalOptim.Display = 'iter';
% SMopts.optsGlobalOptim.PopulationSize = 150;
% SMopts.optsGlobalOptim.Generations = 20;

SMopts.globalSolver = GlobalSearch('Display','iter','NumStageOnePoints',1000,'NumTrialPoints',1300);
SMopts.optsGlobalOptim.Display = 'iter';

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
% SMopts.optsLocalOptim.MaxFunEvals = 11;

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

plotOpts = {};
plotOpts.plotModelOpts = {};
plotOpts.plotModelOpts.ylim = [-70 0];
plotOpts.plotModelOpts.pbaspectVec = [1 1 1];

%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit, Sinit, SMopts, Mf, Mc, OPTopts, plotOpts);

keyboard;

