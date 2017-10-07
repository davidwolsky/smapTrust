% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = 12e9;
fmax = 18e9;  
Nm = 41;	% Number of frequencies

% Initial input parameters
%       {'L1',   'L2'}
xinit = [15,     15]';  
% Initial implicit parameters
%       {'eps_r1',  'eps_r2',   'h1',   'h2'}
xpinit = []';

filename = mfilename('SM_GWhybrid_CST_AWR.m');
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '');

% Set up fine model
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'CST\'];
Mf.name = 'GWhybrid_Fine';
Mf.solver = 'CST';
Mf.params = {'L1',  'L2'};
Mf.ximin =  [14,    14]';
Mf.ximax =  [18,    18]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (AWR)
Mc.path = [currentPath,'AWR\'];
Mc.name = 'GWhybrid_Coarse';
Mc.solver = 'AWR';
Mc.params = {'L1',  'L2'};
Mc.ximin = [12,     12]';
Mc.ximax = [18,     18]';
Mc.Iparams ={};
Mc.xpmin =  []';
Mc.xpmax =  []';
Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

%% Set up the SM
% The initial SM structure
Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA = 0;
SMopts.getB = 0;
SMopts.getc = 1;
SMopts.getG = 0;
SMopts.getxp = 0;
SMopts.getF = 0;
SMopts.getE = 0;
SMopts.getd = 0;


%% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 5;
% OPTopts.TRNi = OPTopts.Ni*2;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S1,1','S3,1','S4,1'};
OPTopts.globOpt = 0;
OPTopts.globOptSM = 0;
% OPTopts.globOptSM = 1;
%
OPTopts.goalType =      {'lt','eq','eq'};
OPTopts.goalResType =   {'S1,1_dB', 'S3,1_dB',  'S4,1_dB'};
OPTopts.goalVal =       {-15,       -3,         -3};
OPTopts.goalWeight =    {1,         1,          1};
OPTopts.goalStart =     {14.6e9,    14.8e9,     14.8e9};
OPTopts.goalStop =      {15.4e9,    15.2e9,     15.2e9};
OPTopts.errNorm =       {1,1,1};
% OPTopts.goalType =      {'lt',      'gt',       'lt'};
% OPTopts.goalResType =   {'S2,1_dB', 'S2,1_dB',  'S2,1_dB'};
% OPTopts.goalVal =       {-20,       -3,         -20};
% OPTopts.goalWeight =    {1,         1,          1};
% OPTopts.goalStart =     {4.5e9,     4.9e9,      5.3e9};
% OPTopts.goalStop =      {4.7e9,     5.1e9,      5.5e9};
% OPTopts.errNorm =       {1,         1,          1};

OPTopts.plotIter = 1;
OPTopts.TolX = 10e-2;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
OPTopts.DeltaInit = 0.25;
OPTopts.startWithIterationZero = true;
% OPTopts.prepopulatedSpaceFile = 'SMLog_bandpassFilter.mat'

OPTopts.optsLocalOptim.Display = 'iter-detailed';


%% Set up the parameter extraction
SMopts.globalSolver = 'ga';
SMopts.optsGlobalOptim = optimoptions('ga');
SMopts.optsGlobalOptim.Display = 'final';

SMopts.localSolver = 'fmincon';
SMopts.optsLocalOptim = optimoptions('fmincon');
SMopts.optsLocalOptim.Display = 'iter-detailed';
SMopts.optsLocalOptim.Diagnostics = 'on';
% SMopts.normaliseAlignmentParameters = 0;
% SMopts.optsLocalOptim.DiffMinChange = 1e-5;
SMopts.normaliseAlignmentParameters = 1;
SMopts.optsLocalOptim.DiffMinChange = 1e-4;
SMopts.plotAlignmentFlag = 0;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.optsFminS = optimset('display','iter');
SMopts.optsPBIL.display =  'iter';
SMopts.optsPBIL.Nfeval = 5000;
SMopts.errNorm = 1;
errW = zeros(size(Mf.freq));
errW(Mf.freq > OPTopts.goalStart{1} & Mf.freq < OPTopts.goalStop{1}) = 1;
% errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 1
% errW(Mf.freq > OPTopts.goalStart{3} & Mf.freq < OPTopts.goalStop{3}) = 1
SMopts.errW = errW;
% SMopts.wk = 5;
SMopts.wk = []; % <---
% SMopts.wk = 1;


%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;
