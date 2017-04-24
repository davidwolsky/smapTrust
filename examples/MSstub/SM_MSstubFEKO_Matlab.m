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
xinit = [70e-3]';  % Initial input parameters (only ls in this case)
xpinit = [2.1]';   % Initial implicit parameters (only eps_r in this case)


% Set up fine model
%setenv('PATH', [getenv('PATH') ';C:\Program Files\Altair\14.0\feko\bin']);
Mf.path = 'C:\Users\19718330\Documents\GitHub\masters\ws_smap\smapTrust\examples\MSstub\FEKO\';
%Mf.path = '/home/rib/Documents/masters/smap/smap/examples/MSstub/FEKO/';
Mf.name = 'MSstubOpen';
Mf.solver = 'FEKO';
Mf.params = {'ls'};
% Mf.ximin = [0]';
% Mf.ximax = [1]';
Mf.ximin = [50e-3]';
Mf.ximax = [90e-3]';
Mf.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (MATLAB)
Mc.path = [pwd,'\'];
Mc.path = 'C:\Users\19718330\Documents\GitHub\masters\ws_smap\smapTrust\examples\MSstub\AWR\';
Mc.name = @MSstubCoarse;  % Must pass a function handle if MATLAB is the simulator...
Mc.solver = 'MATLAB';
Mc.params = {'x';'xp';'f'}; % Must be this order (and names) for MATLAB sims.  Can omit xp or f though...
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.xpmin = [1.8]';
Mc.xpmax = [3.0]';
% Mc.xpmin = [2]';
% Mc.xpmax = [2.2]';
Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);

% Set up coarse model (FEKO)
% Mc.path = 'c:\Users\ddv\Dropbox\Work\MATLAB\SpaceMapping\examples\MSstub\FEKO\';
% Mc.path = 'C:\Users\19718330\Documents\GitHub\masters\ws_smap\smapTrust\examples\MSstub\FEKO\';
% Mc.name = 'MSstubOpenCoarse';
% Mc.solver = 'FEKO';
% Mc.params = {'ls'}; 
% Mc.ximin = Mf.ximin;
% Mc.ximax = Mf.ximax;
% Mc.Iparams = {'eps_r'};
% Mc.xpmin = [2]';
% Mc.xpmax = [2.2]';
% Mc.freq = reshape(linspace(fmin,fmax,Nm),Nm,1);


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

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.optsFminS = optimset('display','iter');
SMopts.optsPBIL.display =  'iter';
SMopts.optsPBIL.Nfeval = 5000;
SMopts.errNorm = 1;
% SMopts.wk = 5;
SMopts.wk = 0;

% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 3;
% OPTopts.TRNi = OPTopts.Ni*2;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'S11'};
OPTopts.globOpt = 0;
OPTopts.globOptSM = 1;
OPTopts.goalType = {'minimax'};
% OPTopts.goalType = {'lt'};
% OPTopts.goalResType = {'S11_complex'};
OPTopts.goalResType = {'S11_dB'};
% OPTopts.goalVal = {0.1};
OPTopts.goalVal = {-20};
OPTopts.goalWeight = {1};
OPTopts.goalStart = {1.30e9};
OPTopts.goalStop = {1.45e9};
OPTopts.errNorm = {1};
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

%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;

