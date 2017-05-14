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
% Initial input parameters
%       {'L1',  'L2',   'L3',   'L4',   'g'}
xinit = [6.784, 4.890,  6.256,  5.280,  0.0956]';  
% Initial implicit parameters
% TODO_DWW: follow through 
xpinit = [2.1]';

keyboard
currentPath = pwd

% Set up fine model
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = [currentPath,'\FEKO\'];
Mf.name = 'BandpassFilter';
Mf.solver = 'FEKO';
Mf.params = {'L1',  'L2',   'L3',   'L4',   'g'};
Mf.ximin =  [6.0,   4.0,    6.0,    4.5,    0.085]';
Mf.ximax =  [7.0,   5.0,    7.0,    5.5,    0.100]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (AWR)
Mc.path = [currentPath,'\AWR\'];
Mc.name = 'BandpassFilter';
Mc.solver = 'AWR';
Mc.params = {'L1',  'L2',   'L3',   'L4',   'g'};
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.Iparams ={'eps_r1',  'eps_r2',   'h1',   'h2'};
Mc.xpmin =  [8.5,       8.5,        0.60,   0.60]';
Mc.xpmax =  [9.5,       9.5,        0.70,   0.70]';
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
OPTopts.SParamMaxPortNumber = 3;
OPTopts.globOpt = 0;
OPTopts.globOptSM = 1;
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
OPTopts.goalStart = {1.40e9, 1.60e9};
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
OPTopts.TolX = 10e-2;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
OPTopts.DeltaInit = 0.25;
OPTopts.testEnabled = true;

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
SMopts.errW = errW;
% SMopts.wk = 5;
SMopts.wk = 0;


%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;
