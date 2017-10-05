function [Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts)

% Space Mapping main loop

% Inputs:
% xinit:    Initial values of input parameters [Nn,1]
% Sinit:    Initial SM structure (see buildSurr.m/evalSurr.m)
%           Important to include initial implicit parameters
% SMopts:   SM options (see buildSurr.m/evalSurr.m)
% Mf:       Fine model structure
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'FEKO'/'MATLAB' (for now)
%   params:     Cell array of parameter names - same order as xinit {Nn,1}
%   The following (2) limits are only for warning generation - not used in optimization
%   ximin:  Vector of minimum values to limit the parameters [Nn,1] 
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% Mc:       Coarse model structure (can be cell array of structures if more than one type has to be calculated to get all the fine model responses)
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'FEKO'/'MATLAB' (for now)
%   params:     Cell array of parameter names - same order as xinit {Nn,1}
%   The following (2) limits are only for warning generation - not used in optimization
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   Iparams:     Cell array of implicit parameter names - same order as xinit {Nn,1}
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% OPTopts:  Optimization options structure for external loop
%   ximin:      Vector of minimum values to limit the parameters [Nn,1] (Compulsory)
%   ximax:      Vector of maximum values to limit the parameters [Nn,1] (Compulsory)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types (so far):
%               'Sb,a' - s-parameters are handled internally using both real and imaginary parts.
%               'Gen' - generic case for use with MATLAB models.
%   Ni:         Maximum number of iterations
%   TRNi:       Maximum number of iterations for the Trust region loop.
%               To turn the TR off use TRNi=1.
%   globOpt:    Flag to run a globabl optimisation routine (1 for only first iteration, 
%               2 for all iterations) (default 0)
%   globOptSM:  Flag to run the global optimisation routine during the PE process 
%               (1 for only first iteration, 2 for all iterations) (default 0)
%   globalSolver:       String to choose the global optimiser solver type (default 'ga').
%   optsGlobalOptim:    Problem options for the global optimiser chosen (e.g. optimoptions('ga')).
%                       This needs to be specified when a non-default solver is chosen.
%   localSolver:        String to choose the local optimiser solver type (default 'fmincon').
%   optsLocalOptim:     Problem options for the local optimiser, e.g. optimoptions('fmincon'). 
%                       This needs to be specified when a non-default solver is chosen.
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'Sb,a_complex'
%               'Sb,a_real'
%               'Sb,a_imag'
%               'Sb,a_dB'
%               'Sb,a_abs'
%               'Sb,a_angle'
%               'Sb,a_deg' - degrees
%               'Gen' - general            
%   goalType:   Cell array of goal types {1,Ng}
%               Valid types:
%                   'lt' (Less than)
%                   'gt' (Greater than)
%                   'eq' 
%                   'minimax' 
%                   'bw' (Like 'lt' but also maximizes bandwidth - requires goalCent value)
%   goalVal:    Cell array of goal values {1,Ng} - same order as goalType
%   goalWeight: Vector of goal weights [1,Ng]
%   goalStart:  Cell array of start of valid goal domain {1,Ng}
%   goalStop:   Cell array of stop of valid goal domain {1,Ng}
%   goalCent:   Cell array of centre point of goal domain {1,Ng} (used by the 'bw' goalType) (optional)
%   errNorm:    Cell array of type of error norms to use for optimization {1,Ng}
%               Valid types: (1,2,inf)
%   TolX:       Termination tolerance on X [ positive scalar - default 10^-2]
%   plotIter:   Flag to plot the responses after each iteration
%   eta1:       A factor used by the TR to define the bound governing when to keep or reduce the radius.
%   eta2:       A factor used by the TR to decide the bound governing when to keep or grow the radius.
%   alp1:       A factor used by the TR to define the rate at which the radius grows with a very successful run.
%   alp2:       A factor used by the TR to define the rate at which the radius shrinks for divergent fine and surrogate runs.
%   DeltaInit:  The initial trust region radius.
%   startWithIterationZero: Starts with the xinit value right away. This is typrically used for reproducing test cases
%                           where it is not desired for an initial optimisation phase to take place.  
%   prepopulatedSpaceFile: File name used to preload fine model runs so that a better surrogate can be calculated.

% Returns:
% Ri:   Structure containing the responses at each iteration
% Si:   Structure containing the surrogates at each iteration
% Pi:   Structure containing the parameters at each iteration
% Ci:   Structure containing the costs at each iteration
% Oi:   Structure containing the 'optimizer' information at each iteration
% Li:   Structure containing the limiting information at each iteration
% Ti:   Structure containing the trust region information for each iteration and fine model evaluations 

% Date created: 2015-03-06
% Dirk de Villiers and Ryno Beyers
% Last Modified: 2016-10-17
% Updates:
% 2015-03-06: Write function shell and basic functionality
% 2015-03-09: Continue with shell and basic functionality
% 2015-03-10: Continue with shell and basic functionality
% 2015-03-11: Continue with shell and basic functionality
% 2015-03-12: Continue with shell and basic functionality
% 2015-03-13: Continue with shell and basic functionality
% 2015-03-14: Add MATLAB case to fine model evaluation function
% 2015-03-27: Add costeq cost function
%             Add plotIter functionality
% 2015-03-30: Add FEKO functionality to the fine model function
% 2015-04-04: Add FEKO functionality in the coarse model function
% 2015-05-07: Add limits for x/xp in the cost function
% 2015-06-26: Update fineMod to rather switch between solvers
%             Update coarseMod to rather switch between solvers 
%             Removed costFunc and dependants from function to be used
%             outside as well
% 2016-02-24: Fix some issues with xp limits if user supplies them and does
%             not want xp SM included
% 2016-03-13: Update how search space box limits are handled - they are
%             compulsory
% 2016-03-14: Add the delta_x exit criterion
% 2016-05-27: Normalize data for optimization
% 2016-08-10: Normalization fixed because it was terrible!
% 2016-08-21: Re-factored the main loop to get fine model evaluation at the end and first 
%             iteration setup before loop.
% 2016-10-17: Introduced the basic trust region (BTR) based on Trust-Region Methods by 
%             A. R. Conn, N. I. M. Gould and P. L. Toint   
% 2016-10-17: Create a test suite that is relatively deterministic for a mathematical model. 
% 2016-11-13: Allow the TR to be turned off using TRNi=1.
% 2017-02-13: Bug fix and restructuring related to correcting the surrogate model being used 
%             at successive TR iterations.
% 2017-03-01: Introduced a TRterminate flag to ensure that an unsuccessful run is not plotted 
%             and that counts are not incremented incorrectly in this case.
% 2017-06-27: Log files now capture the core variable workspace after each successful iteration.
%             This space can then be read in before the start of the main loop. So far only the 
%             Ti.xi_all, Ti.Rfi_all, Ti.costF_all are loaded from the previous space. They are 
%             used to build a better surrogate. In future a surrogate could potentially be 
%             specified from outside SMmain. The 'prepopulatedSpaceFile' variable on OPTopts 
%             is used to specify the space that should be loaded. 
% 2017-08-14: Added validation for the coarse AWR model evaluation to pick up errors in simulations.
%             Also fixed coarse models being clamped to the values that are acceptable in the evaluation.
%             The optimisation can give output beyond what is defined as the constraints.

% ===== ToDo =====
% - Global optimisor
% - Next example
% - Check ToDo and CRC_ in the code
% =====      =====

% Set defaults
Ni = 10;    % Maximum number of iterations
TRNi = Ni;  % Maximum number of iterations for the Trust region loop
TolX = 10^-2;
TolXnorm = 0;
globOpt = 0;
localSolver = 'fmincon';
optsLocalOptim = optimoptions('fmincon');
globalSolver = 'ga';
optsGlobalOptim = optimoptions('ga');
globOptSM = 0;
plotIter = 1;
eta1 = 0.05;
eta2 = 0.9;
alp1 = 2.5;
alp2 = 0.25;
startWithIterationZero = 0;
prepopulatedSpaceFile = '';
DeltaInit = 0.25;

if isfield(OPTopts,'Ni'), Ni = OPTopts.Ni; end
if isfield(OPTopts,'TRNi'), TRNi = OPTopts.TRNi; end
if isfield(OPTopts,'TolX'), TolX = abs(OPTopts.TolX); end % Force positive
if isfield(OPTopts,'globOpt'), globOpt = OPTopts.globOpt; end
if isfield(OPTopts,'globOptSM'), globOptSM = OPTopts.globOptSM; end
if isfield(OPTopts,'localSolver'), localSolver = OPTopts.localSolver; end
if isfield(OPTopts,'optsLocalOptim'), optsLocalOptim = OPTopts.optsLocalOptim; end
if isfield(OPTopts,'globalSolver'), globalSolver = OPTopts.globalSolver; end
if isfield(OPTopts,'optsGlobalOptim'), optsGlobalOptim = OPTopts.optsGlobalOptim; end
if isfield(OPTopts,'plotIter'), plotIter = OPTopts.plotIter; end
if isfield(OPTopts,'eta1'), eta1 = OPTopts.eta1; end
if isfield(OPTopts,'eta2'), eta2 = OPTopts.eta2; end
if isfield(OPTopts,'alp1'), alp1 = OPTopts.alp1; end
if isfield(OPTopts,'alp2'), alp2 = OPTopts.alp2; end
if isfield(OPTopts,'DeltaInit'), DeltaInit = OPTopts.DeltaInit; end
if isfield(OPTopts,'startWithIterationZero'), startWithIterationZero = OPTopts.startWithIterationZero; end
if isfield(OPTopts,'prepopulatedSpaceFile'), prepopulatedSpaceFile = OPTopts.prepopulatedSpaceFile; end

% Set up models - bookkeeping
Nq = 0;
Nn = length(xinit);
xi{1} = reshape(xinit,Nn,1);
if isfield(Sinit,'xp')
    Nq = length(Sinit.xp);
    xpi{1} = reshape(Sinit.xp,Nq,1);
else
    [Sinit.xp] = deal([]);
end
Nr = length(OPTopts.Rtype); % Number of responses requested
Mc.Rtype = OPTopts.Rtype;
Mf.Rtype = OPTopts.Rtype;
Sinit.M = Mc;
Sinit.coarse = @coarseMod;
if isfield(Mc,'freq')
    Sinit.f = Mc.freq;
    fc = Mc.freq;
else
    fc = [];
end
useAllFine = 0;
if isfield(SMopts,'wk')
    % Assume user wants to use all calculated fine responses to fit the surrogate...
    useAllFine = 1;
end


% !!! NOT SURE IF THIS PART SHOULD BE IN - WORKS FINE WITHOUT... 
% 
% % Limit the parameter space
% if all([isfield(Mc,'ximin'),isfield(Mc,'ximax')])
%     Sinit.ximin = Mc.ximin;
%     Sinit.ximax = Mc.ximax;
% else
%     error('ximin and ximax required for the coarse model in Mc');
% end
% if isfield(Sinit,'xp')
%     if isfield(SMopts,'getxp') 
%         if SMopts.getxp 
%             if all([isfield(Mc,'xpmin'),isfield(Mc,'xpmax')])
%                 Sinit.xpmin = Mc.xpmin;
%                 Sinit.xpmax = Mc.xpmax;
%             else
%                 error('xpmin and xpmax required for the coarse model in Mc');
%             end
%         end
%     elseif isfield(SMopts,'getG')
%         if SMopts.getG > 0 
%             if all([isfield(Mc,'xpmin'),isfield(Mc,'xpmax')])
%                 Sinit.xpmin = Mc.xpmin;
%                 Sinit.xpmax = Mc.xpmax;
%             else
%                 error('xpmin and xpmax required for the coarse model in Mc');
%             end
%         end
%     end
% end


% Enter the main loop
%   0)  Normalise parameters
%       * Optimise coarse model to find initial alignment position
%       * Evaluate the fine model at starting position
%       * Get the initial response
%   ->  Test for convergence, main loop count limit and tolerance:
%       1)  Use TR to set up bounds for the optimiser.
%       ->  Check for TR iteration success, TR count limit and tolerance:
%           1)  Optimize the current model Si{ii} to find the next test point xi (xi{ii+1}).
%               Evaluation are placed forward into the next (ii+1) iteration place-holder. 
%               This is over-written on the next iteration if the runs is not successful 
%               i.e. the surrogate and fine models diverge. 
%           2)  Evaluate the fine model at next position (Rfi{ii+1}), this will be overwritten if 
%               the TR iteration is not successful.
%           3)  Calculate the costs change for the fine model and the surrogates (between ii and ii+1). 
%               The difference is stored in rho{ii} and clipped to zero if one, or both, costs get worse.
%           4)  Calculate the step (sk{ii}) between the current and next parameter.
%           5)  Decide if this step is successful.
%           6)  Depending on how successful either keep the current normalised radius (Deltan{ii}) or grow it. 
%               If it is unsuccessful (the costs diverge) then shrink the radius and try this step again 
%               (increase kk but not ii).
%           7)  If this run was successful set the target count to the next position ii+1, else leave it at ii.
%               This is done to ensure that the surrogate, that is going to be built with the new fine model
%               evaluation, is placed in the correct place. If the run was not successful then the ii 
%               position gets updated and the optimisation has more data to use for its next evaluation.
%           ->  Iterate over the responses: 
%               *   Get the response Rsi{targetCount} of the current iteration surrogate (Si{ii}) at the new
%                   test point (xi{ii+1}). 
%               *   Re-evaluate the current surrogate (Si{ii} with xi{ii}) and that for the next step (Si{ii+1} 
%                   with xi{ii+1}). This is done with the new fine model evaluation included. The re-evaluation 
%                   is also done so that the surrogate costs can be compared correctly.
%               *   Align the model (Rsai) at using the new surrogate at the next point (xi{ii+1}).

specF = 0;  % Flag to test if the fine model reached spec
TolX_achieved = 0;
TRterminate = 0;
space = {};
[limMin_f{1},limMax_f{1},limMin_c{1},limMax_c{1}] = deal(zeros(size(xi{1})));

% Set to zero before initialisation depending on starting set
ii = 0;
% Normalize the optimization parameters
ximinn = OPTopts.ximin - OPTopts.ximin;
ximaxn = OPTopts.ximax./OPTopts.ximax;
% xinitn = (xinit - OPTopts.ximin)./(OPTopts.ximax - OPTopts.ximin);

% The initial trust region radius
Ti.Deltan{1} = DeltaInit;
Ti.Delta{1} = DeltaInit.*(OPTopts.ximax - OPTopts.ximin);

% For the initial starting point ii=1
ii = 1;
% Normalize the optimization parameters
ximinn = OPTopts.ximin - OPTopts.ximin;
ximaxn = OPTopts.ximax./OPTopts.ximax;
% CRC_DDV: DWW: We need to come up with better names for these xinitn, init, xin.
% camel case or paskal or with underscores or something >.< 
xinitn = (xinit - OPTopts.ximin)./(OPTopts.ximax - OPTopts.ximin);
% The initial trust region radius
Ti.Deltan{1} = DeltaInit;
Ti.Delta{1} = DeltaInit.*(OPTopts.ximax - OPTopts.ximin);

if ~startWithIterationZero
    % Optimize coarse model to find initial alignment position
    problem = {};
    problem.x0 = xinitn;
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = ximinn;
    problem.ub = ximaxn;
    problem.nonlcon = [];
    if globOpt
        problem.fitnessfcn = @(tempXinGlobal) costSurr(tempXinGlobal,Sinit,OPTopts)
        problem.nvars = length(problem.x0);
        problem.solver = globalSolver;
        problem.options = optsGlobalOptim;
        [xinGlobal,fval,exitflag,output] = doOptimisation(problem);
        xinGlobal = reshape(xinGlobal, length(xinGlobal),1)
        % Start with global search to get initial value.
        problem.x0 = xinGlobal;
    end
    problem.objective = @(tempXin) costSurr(tempXin,Sinit,OPTopts)
    problem.solver = localSolver;
    problem.options = optsLocalOptim;
    [xin{1}, costS{1}, exitflag, output] = doOptimisation(problem)
else
    startWithIterationZero
    xin{1} = xinitn;
end
% De-normalize input vector
xi{1} = xin{1}.*(OPTopts.ximax - OPTopts.ximin) + OPTopts.ximin;

% setup initialisation set/space
if length(prepopulatedSpaceFile) > 1
    space = getPreppopulatedSpaceFrom(prepopulatedSpaceFile);
    assert(size(space.Rfi{1}{1}.f,2) == Nr, 'The number of results expected and those imported must correspond.')
    assert(size(space.Ti.xi_all,1) == size(space.Ti.Rfi_all,1), 'The number of fine model runs and parameters must correspond')

    Ti.xi_all = space.Ti.xi_all;
    Ti.Rfi_all = space.Ti.Rfi_all;
    Ti.costF_all = space.Ti.costF_all;

    Rci{1} = coarseMod(Mc, xi{1}, Sinit.xp, fc);
    Rfi{1} = fineMod(Mf, xi{1});

    Ti.xi_all{end+1} = xi{1};
    Ti.Rfi_all{end+1} = Rfi{1};

    for rr = 1:Nr
        r = {};
        for iii = 1:length(Ti.Rfi_all)
            r{end+1} = Ti.Rfi_all{iii}{rr}.r;
        end
        if globOptSM > 0, SMopts.globOpt = 1; end
        Si{1}{rr} = buildSurr(Ti.xi_all,r,Sinit,SMopts);
        Rsi{1}{rr}.r = evalSurr(xi{1},Si{1}{rr});
        Rsi{1}{rr}.t = Rci{1}{rr}.t;
        if isfield(Rci{1}{rr},'f'), Rsi{1}{rr}.f = Rci{1}{rr}.f; end
        Rsai{1}{rr} = Rsi{1}{rr};
    end
    costS{1} = costSurr(xin{1},Si{1}{:},OPTopts);

    Ti.successCount = [1];
    Ti.Si_all{1} = Si{1};
    Ti.costS_all{1} = costS{1};
    Ti.costChangeS{1} = costS{1};
    Ti.rho = [];
    Ti.rho_all = [];

    costF{1} = costFunc(Ti.Rfi_all{end},OPTopts);
    Ti.costF_all{end+1} = costF{1};
    Ti.costChangeF{1} = costF{1}; 

else 
    % TODO_DWW: Leave this in and formalise... 
    display(['--- TODO_DWW: Initialising ---'])
    Rci{1} = coarseMod(Mc,xi{1},Sinit.xp,fc);
    Rfi{1} = fineMod(Mf,xi{1});
    for rr = 1:Nr
        if globOptSM > 0, SMopts.globOpt = 1; end
        Si{1}{rr} = buildSurr(xi{1},Rfi{1}{rr}.r,Sinit,SMopts);
        Rsi{1}{rr}.r = evalSurr(xi{1},Si{1}{rr});
        Rsi{1}{rr}.t = Rci{1}{rr}.t;
        if isfield(Rci{1}{rr},'f'), Rsi{1}{rr}.f = Rci{1}{rr}.f; end
        Rsai{1}{rr} = Rsi{1}{rr};
    end
    costS{1} = costSurr(xin{1},Si{1}{:},OPTopts);

    Ti.xi_all{1} = xi{1};
    Ti.Rfi_all{1} = Rfi{1};
    Ti.successCount = [1];
    Ti.Si_all{1} = Si{1};
    Ti.costS_all{1} = costS{1};
    Ti.costChangeS{1} = costS{1};
    Ti.rho = [];
    Ti.rho_all = [];

    costF{1} = costFunc(Rfi{1},OPTopts);
    Ti.costF_all{1} = costF{1};
    Ti.costChangeF{1} = costF{1}; 
end

% Plot the initial fine, coarse, optimised surrogate and aligned surrogate
plotModels(plotIter, ii, Rci, Rfi, Rsi, Rsai, OPTopts);

display(['--- TODO_DWW: Start loop ---'])

while ii <= Ni && ~specF && ~TolX_achieved && ~TRterminate
    % Coming into this iteration as ii now with the fine model run here already and responses available. 
    % Exit if spec is reached (will typically not work for eq and never for minimax, and bw is explicitly excluded)
    % CRC_DDV: DWW: Think there should be some basic specF calculations. 
    % No use different goal less than 
    if costF{ii} == 0 && isempty(find(ismember(OPTopts.goalType,'bw'),1))
        specF = 1;
    else
        specF = 0;
        TRsuccess = 0;
        % TR loop count
        kk = 1; 
        while ~TRsuccess && kk <= TRNi && ~TolX_achieved && ~TRterminate
            % Set up TR boundaries or remove them if TR is not to be used.
            if ( TRNi == 1 )
                ximinnTR = ximinn;
                ximaxnTR = ximaxn;
            else
                ximinnTR = max((xin{ii} - Ti.Deltan{ii}),ximinn);
                ximaxnTR = min((xin{ii} + Ti.Deltan{ii}),ximaxn);
            end
            
            display(['--- TODO_DWW:loop optimisation ', num2str(ii), ' ---'])
            problem = {};
            problem.x0 = xin{ii};
            problem.Aineq = [];
            problem.bineq = [];
            problem.Aeq = [];
            problem.beq = [];
            problem.lb = ximinnTR;
            problem.ub = ximaxnTR;
            problem.nonlcon = [];
            if globOpt == 2
                problem.fitnessfcn = @(tempXinGlobal) costSurr(tempXinGlobal, Si{ii}{:}, OPTopts)
                problem.nvars = length(problem.x0);
                problem.solver = globalSolver;
                problem.options = optsGlobalOptim;
                [xinGlobal,fval,exitflag,output] = doOptimisation(problem);
                xinGlobal = reshape(xinGlobal, length(xinGlobal),1)
                % Start with global search to get initial value.
                problem.x0 = xinGlobal;
            end
            problem.objective = @(tempXin) costSurr(tempXin, Si{ii}{:}, OPTopts)
            problem.solver = localSolver;
            problem.options = optsLocalOptim;
            [xin{ii+1}, costSi, exitflag, output] = doOptimisation(problem)
            assert( costS{ii} >= costSi ); % The cost must have gotten better else chaos!
			Ti.costS_all{end+1} = costSi;

            % De-normalize input vector. The new input vector that is.
            xi{ii+1} = xin{ii+1}.*(OPTopts.ximax - OPTopts.ximin) + OPTopts.ximin;
            Ti.xi_all{end+1}  = xi{ii+1};
            
            % L2 norm describing the parameter space distance between the points
            TolXnorm = norm((xin{ii+1} - xin{ii}),2);
            TolX_achieved = TolXnorm < TolX;
            
            enforceFineModelLimits();
            if ( plotIter )
                Rci{ii+1} = coarseMod(Mc,xi{ii+1},Sinit.xp,fc);
            end

            Rfi{ii+1} = fineMod(Mf,xi{ii+1});
            Ti.Rfi_all{end+1} = Rfi{ii+1};
            
            % Test fine model response
            costF{ii+1} = costFunc(Rfi{ii+1},OPTopts);
            Ti.costF_all{end+1} = costF{ii+1};
            
            % Evaluate results and adjust radius for next iteration
			costChangeF = (costF{ii} - costF{ii+1});
			costChangeS = (costS{ii} - costSi);
            Ti.costChangeF{end+1} = costChangeF; 
            Ti.costChangeS{end+1} = costChangeS; 

			if ( costChangeF > 0 && costChangeS > 0 && abs(costChangeS) > TolX )
				Ti.rho{ii}{kk} = (costChangeF)./(costChangeS);
            else
                Ti.rho{ii}{kk} = 0.0;
            end
            Ti.rho_all(end+1) = Ti.rho{ii}{kk};
			
            Ti.sk{ii} = xin{ii+1}-xin{ii};
            if Ti.rho{ii}{kk} >= eta2
                TRsuccess = 1;
                Ti.Deltan{ii+1} = max(alp1.*norm(Ti.sk{ii}),Ti.Deltan{ii});
                Ti.Delta{ii+1} = Ti.Deltan{ii+1}.*(OPTopts.ximax - OPTopts.ximin);
            elseif Ti.rho{ii}{kk} > eta1
                TRsuccess = 1;
                Ti.Deltan{ii+1} = Ti.Deltan{ii};
                Ti.Delta{ii+1} = Ti.Deltan{ii+1}.*(OPTopts.ximax - OPTopts.ximin);
            else
                TRsuccess = 0;
                Ti.Deltan{ii} = alp2*norm(Ti.sk{ii}); % Shrink current Delta
                Ti.Delta{ii} = Ti.Deltan{ii}.*(OPTopts.ximax - OPTopts.ximin);
                % Keep next iteration clean in-case we are not using the TR
                % or if we encounter a tolerance problem before the next
                % iteration.
                Ti.Deltan{ii+1} = 0;
                Ti.Delta{ii+1} = 0.*(OPTopts.ximax - OPTopts.ximin);
            end

            targetCount = ii;
            if ( TRsuccess || (kk == TRNi) || TolX_achieved )
                targetCount = ii + 1;
            end
            for rr = 1:Nr
                % Get the surrogate response after previous iteration
                % optimization - thus at current iteration position
                Rsi{ii+1}{rr}.r = evalSurr(xi{ii+1},Si{ii}{rr});
                Rsi{ii+1}{rr}.t = Rci{1}{rr}.t;
                if isfield(Rci{1}{rr},'f')
                    Rsi{ii+1}{rr}.f = Rci{1}{rr}.f; 
                end
                if globOptSM < 2, SMopts.globOpt = 0; end
                if ~useAllFine
                    % Re-evaluate the surrogate at the new point. 
                    % CRC_DDV: Not sure about which count to use here...
                    Si{ii+1}{rr} = buildSurr(xi{ii+1},Rfi{ii+1}{rr}.r,Si{ii}{rr},SMopts);
                else
                    % TODO_DWW: Check that it makes sense for this to be inside the loop rr
                    r = {};
                    for iii = 1:length(Ti.Rfi_all)
                        % DISP = ['for - ', num2str(iii),' ',];
                        % disp(DISP)
                        r{end+1} = Ti.Rfi_all{iii}{rr}.r;
                    end
                    % TODO_DWW: decide if we putting this back in
                    % length(Ti.xi_all)
                    % length(r)
                    % assert(length(Ti.xi_all) == length(r), 'The lengths of xi and responses needs to be the same.')
                    
                    display(['--- TODO_DWW: build surr for iteration ', num2str(ii), ' ---'])
                    Si{targetCount}{rr} = buildSurr(Ti.xi_all,r,Si{ii}{rr},SMopts);
                end
                % Also get the currently aligned surrogate for comparison
                Rsai{targetCount}{rr}.r = evalSurr(xi{ii+1},Si{targetCount}{rr});
                Rsai{targetCount}{rr}.t = Rci{1}{rr}.t;
                if isfield(Rci{1}{rr},'f')
                    Rsai{targetCount}{rr}.f = Rci{1}{rr}.f; 
                end
            end % for rr

            Ti.Si_all{end+1} = Si{targetCount};
            costS{targetCount} = costSurr(xin{targetCount},Si{targetCount}{:},OPTopts);
            
            kk = kk+1;  % Increase TR loop count
            
            TRterminate = ~TRsuccess && ( (kk > TRNi) || TolX_achieved );
            % Remove ii+1 entry because it didn't succeed
            if TRterminate
                xi = xi(1:end-1);
                xin = xin(1:end-1);
                Rfi = Rfi(1:end-1);
                costF = costF(1:end-1);
                Rci = Rci(1:end-1);
                Ti.Delta = Ti.Delta(1:end-1);
                Ti.Deltan = Ti.Deltan(1:end-1);

                Rsi = Rsi(1:end-1);
                Rsai = Rsai(1:end-1);
            end

        end % while ~TRsuccess
        
        if (~TRterminate)
            space.ii = ii;
            space.xi = xi;
            space.Rci = Rci;
            space.Rfi = Rfi;
            space.Rsi = Rsi;
            space.Si = Si;
            space.costS = costS;
            space.costF = costF;
            space.limMin_f = limMin_f;
            space.limMax_f = limMax_f;
            space.limMin_c = limMin_c;
            space.limMax_c = limMax_c;
            space.TolX_achieved = TolX_achieved;
            space.Ti = Ti;
            space.TRterminate = TRterminate;

            % TODO_DWW: Remember to update this to a better default
            % createLog('SMLog_MSstub.mat', space);
            createLog('SMLog.mat', space);
%             createLog(['SMLog_',Mf.name,'.mat'], space);
            % Make a (crude) log file
            % save SMLog ii xi Rci Rfi Rsi Si costS costF limMin_f limMax_f limMin_c limMax_c TolX_achieved Ti TRterminate

            % Plot the fine, coarse, optimised surrogate and aligned surrogate
            plotModels(plotIter, ii+1, Rci, Rfi, Rsi, Rsai, OPTopts);
        end
    end
    
    if (~TRterminate)
        ii = ii+1;  % Increase main loop count
        if (ii == Ni)
            display(['====== Terminated due to: max iteration count reached (', num2str(ii), '). ======'])
        end
    else
        % TODO_DWW: refine these messages
        if (specF)
            display(['====== Terminated due to: specF reached. ======'])
        end
        if (TolX_achieved)
            display(['====== Terminated due to: TolX_achieved. ======'])
        end
        if (TRterminate)
            display(['====== Terminated due to: TRterminate. ======'])
        end
    end
end % Main while loop

% Handle output structures
Ri.Rc = Rci;
Ri.Rf = Rfi;
Ri.Rs = Rsi;    % Surrogate after optimization
Ri.Rsa = Rsai;  % Surrogate before optimization, just after alignment at end of previous iteration

Pi = xi;

plotNormalised = true;
plotIterations(true, xi, Ti.Delta, OPTopts, SMopts, Si, plotNormalised, 'Normalised');
plotIterations(true, xi, Ti.Delta, OPTopts, SMopts, Si, ~plotNormalised, 'De-normalised/globalised/universalised');
Ci.costS = costS;
Ci.costF = costF;

Oi.specF = specF;
Oi.TolX_achieved = TolX_achieved;   % Flag
Oi.TolXnorm = TolXnorm; % Actual value
Oi.Ni = ii;
Oi.rho = Ti.rho;
Oi.Delta = Ti.Delta;

Li.limMin_f = limMin_f;
Li.limMax_f = limMax_f;
Li.limMin_c = limMin_c;
Li.limMax_c = limMax_c;

plotCosts(Ti, OPTopts, costS, costF)

% ----- Helper functions -----
function enforceFineModelLimits()
	% Check if fine model is limited 
	if isfield(Mf,'ximin')
		limMin_f{ii} = xi{ii+1} < Mf.ximin;
		if any(limMin_f{ii})
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ', Mf.ximin = ', mat2str(Mf.ximin)) )
		end
	end
	if isfield(Mc,'ximin')
		limMin_c{ii} = xi{ii+1} < Mc.ximin;
		if any(limMin_c{ii})
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ', Mc.ximin = ', mat2str(Mc.ximin)) )
		end
	end
	if isfield(Mf,'ximax')
		limMax_f{ii} = xi{ii+1} > Mf.ximax;
		if any(limMax_f{ii})
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ', Mf.ximax = ', mat2str(Mf.ximax)) )
		end
	end
	if isfield(Mc,'ximax')
		limMax_c{ii} = xi{ii+1} > Mc.ximax;
		if any(limMin_f{ii})
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ',  Mc.ximax = ', mat2str( Mc.ximax)) )
		end
	end
end %enforceFineModelLimits

end % SMmain

% ======================================

function Rf = fineMod(M,xi)

% Rf is a cell array of structures containing the response in Rf.r, the type Rf.t, and the
% (optional) domain (typically frequency) in Rf.f.  Same length as M.Rtype
% xi is an array of input parameters - same order as those specified in M
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB'/'FEKO' (for now)
%               'S1,1_dB'!

% Limit the inputs
if isfield(M,'ximin')
    minI = xi < M.ximin;
    xi(minI) = M.ximin(minI);
    if minI | 0
        warning( strcat('Out of bounds fine model evaluation encountered on ximin = ', ...
            mat2str(M.ximin), ', xi = ', mat2str(xi)) )
    end
    
%   params:     Cell array of parameter names - same order as xinit {Nn,1}
%   ximin:      Vector of minimum values to limit the parameters [Nn,1]
%   ximax:      Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
end
if isfield(M,'ximax')
    maxI = xi > M.ximax;
    xi(maxI) = M.ximax(maxI);
    if maxI | 0
		warning( strcat('Out of bounds fine model evaluation encountered on ximax = ', ...
			mat2str(M.ximax), ', xi = ', mat2str(xi)) )
	end
end

% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end

% Call the correct solver
switch M.solver
    case 'AWR'
        error('AWR solver not implemented yet for fine model evaluations')
    case 'CST'
        Rf = cstMod(M, xi, [], Rtype, []);
    case 'FEKO'
        Rf = fekoMod(M, xi, [], Rtype, []);
    case 'MATLAB'
        Rf = matlabMod(M, xi, [], Rtype, []);
    otherwise
        error(['M.solver unknown for fine model evaluation'])
end

end % fineMod function

% ======================================

function Rc = coarseMod(M, xi, xp, f)

% Rc: is a cell array of structures containing the response in Rc.r, the type Rc.t, and the
% (optional) domain (typically frequency) in Rc.f.  Same length as M.Rtype
% xi: is an array of input parameters - same order as those specified in M
% xp: is an array of implicit parameters - same order as those specified in M
% f: is an array of frequency points where to evaluate the model (optional)
% M: is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'FEKO'/'MATLAB' (for now)
%   params:     Cell array of parameter names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
%               'Sb,a_dB'
%               'Sb,a_complex'
%               'Gen'

% Limit the inputs - this should really never happen...
if isfield(M,'ximin')
    minI = xi < M.ximin;
    if ( any(minI) )
        xi(minI) = M.ximin(minI);
		warning( strcat('Out of bounds coarse model evaluation encountered on ximin = ', ...
			mat2str(M.ximin), ', xi = ', mat2str(xi)) )
    end
end
if isfield(M,'ximax')
    maxI = xi > M.ximax;
    if ( any(maxI) )
        xi(maxI) = M.ximax(maxI);
		warning( strcat('Out of bounds coarse model evaluation encountered on ximax', ...
			mat2str(M.ximax), ', xi = ', mat2str(xi)) )
	end
end
if isfield(M,'xpmin')
    minIp = xp < M.xpmin;
    if ( any(minIp) )
        xp(minIp) = M.xpmin(minIp);
		warning( strcat('Out of bounds coarse model evaluation encountered on xpmin', ...
			mat2str(M.xpmin), ', xi = ', mat2str(xi)) )
    end
end
if isfield(M,'xpmax')
    maxIp = xp > M.xpmax;
    if ( any(maxIp) )
        xp(maxIp) = M.xpmax(maxIp);
		warning( strcat('Out of bounds coarse model evaluation encountered on xpmax', ...
			mat2str(M.xpmax), ', xi = ', mat2str(xi)) )
    end
end

% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end

% Call the correct solver
switch M.solver
    case 'CST'
        error('CST solver not implemented yet for coarse model evaluations')
    case 'FEKO'
        Rc = fekoMod(M, xi, [], Rtype, []);
    case 'AWR'
        Rc = awrMod(M, xi, xp, Rtype, []);
    case 'ADS'
        error('ADS solver not implemented yet for coarse model evaluations')
    case 'MATLAB'
        Rc = matlabMod(M, xi, xp, Rtype, f);
    otherwise
        error(['M.solver unknown for coarse model evaluation'])
end

end % coarseMod function

% ======================================

function cost = costSurr(xin,S,OPTopts)

% Function to run the Surrogate model S at x to return the cost as
% calculated from the information in OPTopts.
% S can be a cell array of Surrogates, each returning a different response
% at xin, with the length of the array equal to the length of the type vector
% in OPTopts.Rtype

if length(S) == 1 && ~iscell(S), S = {S}; end
Nr = length(S); % Number of responses
xn = reshape(xin,length(xin),1);

Rs = cell(1,Nr);
for rr = 1:Nr
    % De-normalize
%     x = xn.*(S{rr}.ximax - S{rr}.ximin) + S{rr}.ximin;
    x = xn.*(OPTopts.ximax - OPTopts.ximin) + OPTopts.ximin;
    Rs{rr}.r = evalSurr(x,S{rr});
    Rs{rr}.t = OPTopts.Rtype{rr};
    if isfield(S{rr},'f'), Rs{rr}.f = S{rr}.f; end
end
cost = costFunc(Rs,OPTopts);

end % costSurr function

% ======================================

function plotCosts(Ti, OPTopts, costS, costF)
markerstr = 'xso+*d^v><ph.xso+*d^v><ph.';
colourstr = 'kbrgmcykbrgmckbrgmcykbrgmc';

figure()
Ni = length(Ti.xi_all);    % Number of iterations
% costSi1 = [];
% for nn = 1:Ni
%     for cc = 1:Ni
%         costSi1(cc) = costSurr(Ti.xin_all{cc},Ti.Si_all{nn}{:},OPTopts);
%     end
%     subplot(2,2,1)
%     plot(costSi1, strcat(markerstr(nn),colourstr(nn)),'LineWidth',2,'MarkerSize',10), grid on, hold on
%     ylabel('costSi')
%     xlabel('Iterations all')
% end
% legend('show') 
% title('costs')
subplot(2,2,[1 2])
if strcmp(OPTopts.goalResType{1},'Gen')
    plot(cell2mat(Ti.costF_all), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'complex')
    plot(real(cell2mat(Ti.costF_all)), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plot(imag(cell2mat(Ti.costF_all)), strcat(markerstr(2),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    legend('real','imag');
    plotRealGoals()
    plotImagGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'dB')
    plot(cell2mat(Ti.costF_all), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
else
    error(['goalResType (', OPTopts.goalResType{1}, ')not found']);
end
ylabel('Ti.costF\_all')
xlabel('Iterations all')
title('Ti.costF\_all')

% Meaningless -> just for debugging
subplot(2,2,3)
if strcmp(OPTopts.goalResType{1},'Gen')
    plot(cell2mat(costS), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'complex')
    plot(real(cell2mat(costS)), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plot(imag(cell2mat(costS)), strcat(markerstr(2),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    legend('real','imag')
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'dB')
    plot(cell2mat(costS), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
else
    error(['goalResType (', OPTopts.goalResType{1}, ')not found']);
end
ylabel('costS')
xlabel('Iterations success points')
title('costS - fairly meaningless')

subplot(2,2,4)
if strcmp(OPTopts.goalResType{1},'Gen')
    plot(cell2mat(costF), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'complex')
    plot(real(cell2mat(costF)), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plot(imag(cell2mat(costF)), strcat(markerstr(2),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    legend('real','imag')
    plotRealGoals()
    plotImagGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'dB')
    plot(cell2mat(costF), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
else
    error(['goalResType (', OPTopts.goalResType{1}, ')not found']);
end
ylabel('costF')
xlabel('Iterations success points')
title('costF')


% ----- goals -----
function plotGoals()
if ( isfield(OPTopts, 'goalVal') )%&& isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        plot([1, Ni], (OPTopts.goalVal{gg})*ones(1,2), 'm', 'LineWidth',2)
    end
end % if validation
end % plotGoals function

function plotRealGoals()
if ( isfield(OPTopts, 'goalVal') )%&& isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        plot([1, Ni], real(OPTopts.goalVal{gg})*ones(1,2), 'm', 'LineWidth',2)
    end
end % if validation
end % plotRealGoals function

function plotImagGoals()
if ( isfield(OPTopts, 'goalVal') )%&& isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        if ( ~isreal(OPTopts.goalVal{gg}) )
            plot([1, Ni], imag(OPTopts.goalVal{gg})*ones(1,2), 'c', 'LineWidth',2)
        end % imag part
    end
end % if validation
end % plotImagGoals function
% -----  -----

end % plotCosts

% ======================================