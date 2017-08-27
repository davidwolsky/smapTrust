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

% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 3;
% OPTopts.TRNi = OPTopts.Ni*2;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S1,1'};
% TODO_DWW: name this nicely and follow through with an assertion
OPTopts.globOpt = 0;
OPTopts.globOptSM = 0;
% OPTopts.globOptSM = 1;
%
% OPTopts.goalType = {'minimax'};
% OPTopts.goalResType = {'S1,1_complex'};
% OPTopts.goalVal = {0.1-1j*0.2};
% OPTopts.goalWeight = {1};
% OPTopts.goalStart = {1.30e9};
% OPTopts.goalStop = {1.45e9};
% OPTopts.errNorm = {1};
%
OPTopts.goalType = {'lt', 'gt'};
OPTopts.goalResType = {'S1,1_dB', 'S1,1_dB'};
OPTopts.goalVal = {-20, -10};
OPTopts.goalWeight = {1, 1};
OPTopts.goalStart = {1.30e9, 1.60e9};
OPTopts.goalStop = {1.45e9, 1.75e9};
OPTopts.errNorm = {1,1};
%
OPTopts.optsPBIL.display =  'iter'; 
OPTopts.optsPBIL.Nfeval = 5000;
OPTopts.optsPBIL.Nbest = 10; % DOM
OPTopts.M_PBIL = 6;
OPTopts.optsFminS = optimset('display','iter');
% OPTopts.optsFminS = optimset('MaxFunEvals',10,'display','iter');
OPTopts.plotIter = 1;
% Normalised tolerance in main loop
OPTopts.TolX = 0.1;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
OPTopts.DeltaInit = 0.25;
OPTopts.testEnabled = true;
% OPTopts.prepopulatedSpaceFile = '.mat';

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.optsFminS = optimset('display','iter');
SMopts.optsPBIL.display =  'iter';
SMopts.optsPBIL.Nfeval = 5000;
% SMopts.errNorm = 1;
SMopts.errNorm = 2;
% errW = 1
errW = zeros(size(Mf.freq));
errW(Mf.freq > OPTopts.goalStart{1} & Mf.freq < OPTopts.goalStop{1}) = 1;
errW(Mf.freq > OPTopts.goalStart{2} & Mf.freq < OPTopts.goalStop{2}) = 1;
SMopts.errW = errW;
% SMopts.wk = 5;
SMopts.wk = 0;
SMopts.plotAlignmentFlag = 1;

%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;
