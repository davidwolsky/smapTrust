% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = 4.5e9;
fmax = 5.5e9;  
Nm = 41;	% Number of frequencies

% Initial input parameters
%       {'L1',  'L2',   'L3',   'L4',   'g'}
xinit = [6.784, 4.890,  6.256,  5.280,  0.0956]';  
% Initial implicit parameters
%       {'eps_r1',  'eps_r2',   'h1',   'h2'}
xpinit = [9.0,      9.0,        0.66,   0.66]';

filename = mfilename('SM_BandpassFilter.m');
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '');

% Set up fine model
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'FEKO\'];
Mf.name = 'bandpassFilter';
Mf.solver = 'FEKO';
Mf.params = {'L1',  'L2',   'L3',   'L4',   'g'};
Mf.ximin =  [6.0,   4.0,    6.0,    4.5,    0.0661]';
Mf.ximax =  [7.0,   5.0,    7.0,    5.5,    0.250]';
% Mf.ximin =  [6.0,   4.0,    6.0,    4.5,    0.085]';
% Mf.ximax =  [7.0,   5.0,    7.0,    5.5,    0.100]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (AWR)
Mc.path = [currentPath,'AWR\'];
Mc.name = 'bandpassFilter';
Mc.solver = 'AWR';
Mc.params = {'L1',  'L2',   'L3',   'L4',   'g'};
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.Iparams ={'eps_r1',  'eps_r2',   'h1',   'h2'};
Mc.xpmin =  [8.0,       8.0,        0.25,   0.50]';
Mc.xpmax =  [10.0,      10.0,       0.661,   1.30]';
% Mc.xpmin =  [6.0,       6.0,        0.50,   0.50]';
% Mc.xpmax =  [10.0,      10.0,       1.30,   1.30]';
% -------
% Microstrip gap EM Quasi-Static: MGAP2
% S = g, W = 0.6, H = h1
% 0.5 < W/H < 2.5
% 0.1 < S/H < 1
% 1 < εr < 216 
% -------
Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up the SM
% The initial SM structure
Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA = 0;
SMopts.getB = 0;
SMopts.getc = 1;
SMopts.getG = 0;
SMopts.getxp = 1;
SMopts.getF = 0;
SMopts.getE = 1;
SMopts.getd = 0;


% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 5;
% OPTopts.TRNi = OPTopts.Ni*2;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S2,1'};
OPTopts.globOpt = 0;
OPTopts.globOptSM = 0;
% OPTopts.globOptSM = 1;
%
OPTopts.goalType =      {'lt',      'gt'};
OPTopts.goalResType =   {'S2,1_dB', 'S2,1_dB'};
OPTopts.goalVal =       {-20,       -3};
OPTopts.goalWeight =    {1,         1};
OPTopts.goalStart =     {4.5e9,     4.9e9};
OPTopts.goalStop =      {4.7e9,     5.1e9};
OPTopts.errNorm =       {1,         1};
% OPTopts.goalType =      {'lt',      'gt',       'lt'};
% OPTopts.goalResType =   {'S2,1_dB', 'S2,1_dB',  'S2,1_dB'};
% OPTopts.goalVal =       {-20,       -3,         -20};
% OPTopts.goalWeight =    {1,         1,          1};
% OPTopts.goalStart =     {4.5e9,     4.9e9,      5.3e9};
% OPTopts.goalStop =      {4.7e9,     5.1e9,      5.5e9};
% OPTopts.errNorm =       {1,         1,          1};
%
OPTopts.plotIter = 1;
OPTopts.TolX = 10e-2;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
OPTopts.DeltaInit = 0.25;
OPTopts.startWithIterationZero = true;
% OPTopts.prepopulatedSpaceFile = 'SMLog_bandpassFilter.mat'

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.errNorm = 1;
% SMopts.errNorm = 2;
% errW = 1;
errW = zeros(size(Mf.freq));
errW(Mf.freq > OPTopts.goalStart{1} & Mf.freq < OPTopts.goalStop{1}) = 1
errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 1
% errW(Mf.freq > OPTopts.goalStart{3} & Mf.freq < OPTopts.goalStop{3}) = 1
SMopts.errW = errW;
% SMopts.wk = 5;
SMopts.wk = 0; % <---
% SMopts.wk = 1;
SMopts.plotAlignmentFlag = 1;


%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;
