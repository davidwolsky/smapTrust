% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = pi/2;
fmax = 7.0*pi;  
Nf = 101;
% f = linspace(0,pi,101);

% Initial input parameters 
xinit = [2, 5.0*pi]';  % Initial input parameters 
xpinit = [2.1]';   % Initial implicit parameters

label = 'SM_MATLABtest_arctan';
filename = mfilename([label, '.m']);
fullpath = mfilename('fullpath');
currentPath = replace(fullpath, filename, '');

% Set up fine model
Mf.path = currentPath;
Mf.label = label;
Mf.name = @fineTest;
Mf.solver = 'MATLAB';
Mf.params = {'x','f'};
Mf.ximin = [1.0,    3.0*pi]';
Mf.ximax = [3.0,    5.0*pi]';
Mf.freq = reshape(linspace(fmin,fmax,Nf),Nf,1);
    
% Set up coarse model (MATLAB)
Mc.path = currentPath;
Mc.label = Mf.label;
Mc.name = @coarseTest;  % Must pass a function handle if MATLAB is the simulator...
Mc.solver = 'MATLAB';
Mc.params = {'x','f'}; % Must be this order (and names) for MATLAB sims.  Can omit xp or f though...
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.xpmin = [2]';
Mc.xpmax = [2.2]';
Mc.freq = reshape(linspace(fmin,fmax,Nf),Nf,1);

% Set up the SM
% The initial SM structure
Sinit = [];
% Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA  = 0;
SMopts.getB  = 0;
SMopts.getc  = 0;
SMopts.getG  = 0;
SMopts.getxp = 0;
SMopts.getF  = 1;
SMopts.getE  = 0;
SMopts.getd  = 0;

% SMopts.globalSolver = 'ga';
% SMopts.optsGlobalOptim = optimoptions('ga');
% SMopts.optsGlobalOptim.Display = 'final';

SMopts.localSolver = 'fmincon';
SMopts.optsLocalOptim = optimoptions('fmincon');
SMopts.optsLocalOptim.Display = 'iter-detailed';
SMopts.optsLocalOptim.Diagnostics = 'on';
% SMopts.optsLocalOptim.DiffMinChange = 1e-6;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
% SMopts.xpmin = Mc.xpmin;
% SMopts.xpmax = Mc.xpmax;
% SMopts.wk = 1.15;
SMopts.plotAlignmentFlag = 1;

% Set up the optimization
OPTopts.startWithIterationZero = 1;
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 3;
OPTopts.TRNi = OPTopts.Ni;
OPTopts.Rtype = {'Gen'};
OPTopts.globOpt = 0;
OPTopts.globOptSM = 0;
OPTopts.onlyDoAlignment = true;

OPTopts.goalType =      {'gt'   };
OPTopts.goalResType =   {'Gen'  };
OPTopts.goalVal =       {1.0    };
OPTopts.goalWeight =    {1.0    };
OPTopts.goalStart =     {10.0    };
OPTopts.goalStop =      {14.0    };
OPTopts.errNorm =       {1      };

% OPTopts.goalType = {'minimax'};
% OPTopts.goalResType = {'Gen'};
% OPTopts.goalVal = {-20};
% OPTopts.goalWeight = {1};
% OPTopts.goalStart = {1.4e9};
% OPTopts.goalStop = {1.6e9};
% OPTopts.errNorm = {1};

%% Run the main loop
[Ri,Si,Pi,Ci,Li] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;
