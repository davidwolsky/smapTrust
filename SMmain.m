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
%               'S11dB'
%               'S11complex'
%               'Gen' - generic case for use with MATLAB models 
%   Ni:         Maximum number of iterations
%   TRNi:       Maximum number of iterations for the Trust region loop.
%               To turn the TR off use TRNi=1.
%   globOpt:    Flag to run PBIL (1 for only first iteration, 2 for all iterations) (default 0)
%   M_PBIL:     Vector of bits for the global search variables (see PBILreal.m)
%   globOptSM:  Flag to run PBIL during the PE process (1 for only first iteration, 2 for all iterations) (default 0)
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'S11dB'
%               'S11complex'
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
%   goalCent:   Cell array of center point of goal domain {1,Ng} (used by the 'bw' goalType) (optional)
%   errNorm:    Cell array of type of error norms to use for optimization {1,Ng}
%               Valid types: (1,2,inf)
%   TolX:       Termination tolerance on X [ positive scalar - default 10^-2]
%   optsFminS:  options for the fminsearch local optimizer
%   optsPBIL:   options for the PBIL global optimizer
%   plotIter:   Flag to plot the responses after each iteration
%   eta1:       A factor used by the TR to define the bound governing when to keep or reduce the radius.
%   eta2:       A factor used by the TR to decide the bound governing when to keep or grow the radius.
%   alp1:       A factor used by the TR to define the rate at which the radius grows with a very successful run.
%   alp2:       A factor used by the TR to define the rate at which the radius shrinks for divergent fine and surrogate runs.
%   DeltaInit:  The initial trust region radius.
% TODO_DWW: Implement this as the next feature for refactoring.  

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
%             Removed costFunc and dependents from function to be used
%             outside as well
% 2016-02-24: Fix some issues with xp limits if user supplies them and does
%             not want xp SM included
% 2016-03-13: Update how search space box limits are handled - they are
%             compulsary
% 2016-03-14: Add the delta_x exit criterion
% 2016-05-27: Normalize data for optimization
% 2016-08-10: Normalization fixed because it was terrible!
% 2016-08-21: Re-factored the main loop to get fine model evaluation at the end and first iteration setup before loop.
% 2016-10-17: Introduced the basic trust region (BTR) based on Trust-Region Methods by A. R. Conn, N. I. M. Gould and P. L. Toint   
% 2016-10-17: Create a test suite that is relatively deterministic for a mathematical model. 
% 2016-11-13: Allow the TR to be turned off using TRNi=1.

% Set defaults
Ni = 10;    % Maximum number of iterations
TRNi = Ni;  % Maximum number of iterations for the Trust region loop
TolX = 10^-2;
globOpt = 0;
M_PBIL = 8;
globOptSM = 0;
optsFminS = optimset('display','none');
optsPBIL = [];
plotIter = 1;
eta1 = 0.05;
eta2 = 0.9;
alp1 = 2.5;
alp2 = 0.25;
testEnabled = 0;
DeltaInit = 0.25;

if isfield(OPTopts,'Ni'), Ni = OPTopts.Ni; end
if isfield(OPTopts,'TRNi'), TRNi = OPTopts.TRNi; end
if isfield(OPTopts,'TolX'), TolX = abs(OPTopts.TolX); end % Force positive
if isfield(OPTopts,'globOpt'), globOpt = OPTopts.globOpt; end
if isfield(OPTopts,'M_PBIL'), M_PBIL = OPTopts.M_PBIL; end
if isfield(OPTopts,'globOptSM'), globOptSM = OPTopts.globOptSM; end
if isfield(OPTopts,'optsFminS'), optsFminS = OPTopts.optsFminS; end
if isfield(OPTopts,'M_PBIL'), M_PBIL = OPTopts.M_PBIL; end
if isfield(OPTopts,'optsPBIL'), optsPBIL = OPTopts.optsPBIL; end
if isfield(OPTopts,'plotIter'), plotIter = OPTopts.plotIter; end
if isfield(OPTopts,'eta1'), eta1 = OPTopts.eta1; end
if isfield(OPTopts,'eta2'), eta2 = OPTopts.eta2; end
if isfield(OPTopts,'alp1'), alp1 = OPTopts.alp1; end
if isfield(OPTopts,'alp2'), alp2 = OPTopts.alp2; end
if isfield(OPTopts,'DeltaInit'), DeltaInit = OPTopts.DeltaInit; end
if isfield(OPTopts,'testEnabled'), testEnabled = OPTopts.testEnabled; end


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
%       - Optimise coarse model to find initial alignment position
%       - Evaluate the fine model at starting position
%       - Get the initial response
%   1)  Test for convergence.
%       2)  Use TR to set up bounds for the optimiser.
%       3)  Optimize the current model Si{ii} to find the next xi (xi{ii+1}).
%       Evaluation are placed forward into the next (ii+1) iteration place-holder. 
%       This is over-written if the runs is not successful, i.e. the surrogate and fine models diverge. 
%       4)  Evaluate the fine model at next position (Rfi{ii+1}).
%       5)  Get the response of the current iteration surrogate (Rsi{ii}) and at the next step (Rsi{ii+1}). 
%       6)  Re-evaluate the current surrogate (Si{ii} with xi{ii}) and that for the next step (Si{ii+1} 
%           with xi{ii+1}). This is done with the new fine model evaluation included. The re-evaluation 
%           is also done so that the surrogate costs can be compared correctly.
%       7)  Align the model at the next step to get Rsai.
%       8)  Calculate the costs change for the fine model and the surrogates (between ii and ii+1). 
%           The difference is stored in rho{ii} and clipped to zero if one, or both, costs get worse.
%       8)  Calculate the step (sk{ii}) between the current and next parameter.
%       9)  Decide if this step is successful.
%      10)  Depending on how successful either keep the current normalised radius (Deltan{ii}) or grow it. 
%           If it is unsuccessful (the costs diverge) then shrink the radius and try this step again 
%           (increase kk but not ii).
%      11)  If successful then clean up the extra fine models that were kept.

specF = 0;  % Flag to test if the fine model reached spec
TolX_achieved = 0;

[limMin_f{1},limMax_f{1},limMin_c{1},limMax_c{1}] = deal(zeros(size(xi{1})));

% For the initial starting point ii=1
ii = 1;
% Normalize the optimization parameters
ximinn = OPTopts.ximin - OPTopts.ximin;
ximaxn = OPTopts.ximax./OPTopts.ximax;
xinitn = (xinit - OPTopts.ximin)./(OPTopts.ximax - OPTopts.ximin);
% The initial trust region radius
Ti.Deltan{1} = DeltaInit;
Ti.Delta{1} = DeltaInit.*(OPTopts.ximax - OPTopts.ximin);

if ~testEnabled
    % Optimize coarse model to find initial alignment position
    if globOpt
        [xinitn,costSi,exitFlag,output] = PBILreal(@(xin) costSurr(xin,Sinit,OPTopts),ximinn,ximaxn,M_PBIL,optsPBIL);
        xinitn = reshape(xinitn,Nn,1);
    end
    LHSmat = [];
    RHSvect = [];
    nonLcon = [];
    [xin{1}, costS{1}] = fminsearchcon(@(xin) costSurr(xin,Sinit,OPTopts),xinitn,ximinn,ximaxn,LHSmat,RHSvect,nonLcon,optsFminS);
else
    testEnabled
    xinitn
    xin{1} = xinitn
    costS{1} = costSurr(xinitn,Sinit,OPTopts)
end


% De-normalize input vector
xi{1} = xin{1}.*(OPTopts.ximax - OPTopts.ximin) + OPTopts.ximin;

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

% TODO: DWW: rename
count_all = 1;
Ti.xi_all{1} = xi{1};
Ti.Rfi_all{1} = Rfi{1};
Ti.successCount = [1];
Ti.Si_all{1} = Si{1};
Ti.costS_all{1} = costS{1};
Ti.rho_all = [];


% Plot the initial fine, coarse, optimised surrogate and aligned surrogate
plotModels(plotIter, 1, Rci, Rfi, Rsi, Rsai, OPTopts);

% Test fine model response
costF{1} = costFunc(Rfi{1},OPTopts);
Ti.costF_all{1} = costF{1};
    
while ii <= Ni && ~specF && ~TolX_achieved
%Coming into this iteration as ii now with the fine model run here already and responses available. 

    % Exit if spec is reached (will typically not work for eq and never for minimax, and bw is explicitly excluded)
    if costF{ii} == 0 && isempty(find(ismember(OPTopts.goalType,'bw'),1))
        specF = 1;
    else
        specF = 0;
        
        % TR
        TRsuccess = 0;
        kk = 1;
        while ~TRsuccess && kk <= TRNi && ~TolX_achieved
            % Set up TR boundaries or remove them if TR is not to be used.
            if ( TRNi == 1 )
                ximinnTR = ximinn;
                ximaxnTR = ximaxn;
            else
                ximinnTR = max((xin{ii} - Ti.Deltan{ii}),ximinn);
                ximaxnTR = min((xin{ii} + Ti.Deltan{ii}),ximaxn);
            end

            % Do optimization
            if globOpt == 2
                [xinitn,costSi,exitFlag,output] = PBILreal(@(xin) costSurr(xin,Si{ii}{:},OPTopts),ximinnTR,ximaxnTR,M_PBIL,optsPBIL);
                xinitn = reshape(xinitn,Nn,1);
            end
            LHSmat = [];
            RHSvect = [];
            nonLcon = [];
            [xin{ii+1}, costSi] = fminsearchcon(@(xin) costSurr(xin,Si{ii}{:},OPTopts),xinitn,ximinnTR,ximaxnTR,LHSmat,RHSvect,nonLcon,optsFminS);
            % De-normalize input vector. The new input vector that is.
            xi{ii+1} = xin{ii+1}.*(OPTopts.ximax - OPTopts.ximin) + OPTopts.ximin;
            
            % L2 norm describing the parameter space distance between the points
            TolXnorm = norm((xin{ii+1} - xin{ii}),2);
            TolX_achieved = TolXnorm < TolX;
            
            enforceFineModelLimits();
            
            if ( plotIter )
                Rci{ii+1} = coarseMod(Mc,xi{ii+1},Sinit.xp,fc);
            end
            Rfi{ii+1} = fineMod(Mf,xi{ii+1});

            Ti.xi_all{end+1}  = xi{ii+1};
			% xin_all{end+1} = xin{ii+1};
            Ti.Rfi_all{end+1} = Rfi{ii+1};
            
            for rr = 1:Nr
                % Get the surrogate response after previous iteration
                % optimization - thus at current iteration position
                % TODO_DWW: comment a bit more and clean up. Explain ii +-1
                Rsi{ii+1}{rr}.r = evalSurr(xi{ii+1},Si{ii+1-1}{rr});
                Rsi{ii+1}{rr}.t = Rci{1}{rr}.t;
                if isfield(Rci{1}{rr},'f')
                    Rsi{ii+1}{rr}.f = Rci{1}{rr}.f; 
                end
                if globOptSM < 2, SMopts.globOpt = 0; end
                if ~useAllFine
                    % Re-evaluate the surrogate at the new point. 
                    Si{ii}{rr}   = buildSurr(xi{ii},Rfi{ii+1}{rr}.r,Si{ii+1-1}{rr},SMopts);
                    Si{ii+1}{rr} = buildSurr(xi{ii+1},Rfi{ii+1}{rr}.r,Si{ii+1-1}{rr},SMopts);
                else
                    % Add prior successful runs
                    iii = 1;
                    r = {};
                    while iii < length(Ti.successCount)
                        DISP = ['while - ', num2str(iii),' ',num2str(Ti.successCount(iii))];
                        disp(DISP)
                        r{iii} = Ti.Rfi_all{Ti.successCount(iii)}{rr}.r;
                        iii = iii+1;
                    end
                    % Additionally add the fine model runs since the last successful run to increase 
                    % data for the surrogate to use.
                    for iii = Ti.successCount(end):length(Ti.Rfi_all)
                        DISP = ['for - ', num2str(iii),' ',];
                        disp(DISP)
                        r{end+1} = Ti.Rfi_all{iii}{rr}.r;
                    end
                    %TODO_DWW: Ti.xi_all should match up with r{}
                    % Re-evaluate the surrogate at the new point. 
                    Si{ii}{rr}   = buildSurr(Ti.xi_all,r,Si{ii+1-1}{rr},SMopts);
                    Si{ii+1}{rr} = buildSurr(Ti.xi_all,r,Si{ii+1-1}{rr},SMopts);
                end
                % Also get the currently aligned surrogate for comparison
                Rsai{ii+1}{rr}.r = evalSurr(xi{ii+1},Si{ii+1}{rr});
                Rsai{ii+1}{rr}.t = Rci{1}{rr}.t;
                if isfield(Rci{1}{rr},'f')
                    Rsai{ii+1}{rr}.f = Rci{1}{rr}.f; 
                end
            end
            Ti.Si_all{end+1} = Si{ii+1};
            
            % Test fine model response
            costF{ii+1} = costFunc(Rfi{ii+1},OPTopts);
			Ti.costF_all{end+1} = costF{ii+1};
			
            % TODO_DWW: Comment here about needing to compare the last and current surrogates
            costS{ii}   = costSurr(xin{ii},Si{ii+1}{:},OPTopts);
            costS{ii+1} = costSurr(xin{ii+1},Si{ii+1}{:},OPTopts);
			Ti.costS_all{end+1} = costS{ii+1};
			
            % Evaluate results and adjust radius for next iteration
			costChangeF = (costF{ii} - costF{ii+1});
			costChangeS = (costS{ii} - costS{ii+1});
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
            
            kk = kk+1;
            count_all = count_all+1;

            % Remove any additional fine model runs and clean up rest of iteration lasting variables.
            if TRsuccess
                Ti.successCount(end+1) = count_all;
                %  TODO_DWW: comment about this.
                %  TODO_DWW: actually, just remove this.
    %             for count = 1:kk-2
    %                 Ti.xi_all(length(Ti.xi_all)-1) = [];
    %                 % xin_all(length(xin_all)-1) = [];
    %                 Ti.Rfi_all(length(Ti.Rfi_all)-1) = [];
    %                 % Ti.Rfi_all(length(Ti.Rfi_all)-1) = [];
    %                 % Si_all(length(Si_all)-1) = [];
    %                 % costS_all(length(costS_all)-1) = [];
    %                 Ti.costF_all(length(Ti.costF_all)-1) = [];
				% end
            end
        end
        
        % Make a (crude) log file
        save SMlog ii xi Rci Rfi Rsi Si costS costF limMin_f limMax_f limMin_c limMax_c Ti
        
        % Plot the fine, coarse, optimised surrogate and aligned surrogate
        plotModels(plotIter, ii+1, Rci, Rfi, Rsi, Rsai, OPTopts);
        
    end
    
    ii = ii+1;
    
end % Main while loop

% Handle output structures
Ri.Rc = Rci;
Ri.Rf = Rfi;
Ri.Rs = Rsi;    % Surrogate after optimization
Ri.Rsa = Rsai;  % Surrogate before optimization, just after alignment at end of previous iteration

Pi = xi;

plotNormalised = true;
plotIterations(true, xi, Ti.Delta, OPTopts, plotNormalised, 'Normalised');
plotIterations(true, xi, Ti.Delta, OPTopts, ~plotNormalised, 'De-normalised/globalised/universalised');
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

% CRC_DWW: DWW: I don't like these function methods either. Can't tell from the 
% function name or statement which iteration we working with, nor do we allow any 
% reuse. Should be passing stuff in and also fn should be out of this loop later. 
function enforceFineModelLimits()
	% Check if fine model is limited 
	if isfield(Mf,'ximin')
		limMin_f{ii} = xi{ii+1} < Mf.ximin;
		if limMin_f{ii} | 0
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ', Mf.ximin = ', mat2str(Mf.ximin)) )
		end
	end
	if isfield(Mc,'ximin')
		limMin_c{ii} = xi{ii+1} < Mc.ximin;
		if limMin_c{ii} | 0
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ', Mc.ximin = ', mat2str(Mc.ximin)) )
		end
	end
	if isfield(Mf,'ximax')
		limMax_f{ii} = xi{ii+1} > Mf.ximax;
		if limMax_f{ii} | 0
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ', Mf.ximax = ', mat2str(Mf.ximax)) )
		end
	end
	if isfield(Mc,'ximax')
		limMax_c{ii} = xi{ii+1} > Mc.ximax;
		if limMin_f{ii} | 0
			warning( strcat('Fine model bound met: xi{ii+1} = ', ...
				mat2str(xi{ii+1}), ',  Mc.ximax = ', mat2str( Mc.ximax)) )
		end
	end
end %enforceFineModelLimits

end % SMmain


function Rf = fineMod(M,xi)

% Rf is a cell array of structures containing the response in Rf.r, the type Rf.t, and the
% (optional) domain (typically frequency) in Rf.f.  Same length as M.Rtype
% xi is an array of input parameters - same order as those specified in M
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB'/'FEKO' (for now)
%               'S11dB' - obvious!

% Limit the inputs
if isfield(M,'ximin')
    minI = xi < M.ximin;
    xi(minI) = M.ximin(minI);
    if minI | 0
        warning( strcat('Out of bounds fine model evaluation encountered on ximin = ', ...
            mat2str(M.ximin), ', xi = ', mat2str(xi)) )
    end
    
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
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

Nn = length(xi);
% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end
Nr = length(Rtype);
Rf = cell(1,Nr);

% Call the correct solver
switch M.solver
    
    case 'CST'
        % Start CST activeX
        cst = actxserver('CSTSTUDIO.Application');
        % Get handle to the model
        mws = invoke(cst,'OpenFile',[M.path,M.name,'.cst']);
        % Update parameters
        for nn = 1:Nn
            invoke(mws,'StoreParameter',M.params{nn},xi(nn));
        end
        invoke(mws,'Rebuild');
        % Run simulation
        solver = invoke(mws,'Solver');
        invoke(solver,'Start');
        
        % Generate output
        for rr = 1:Nr
            if strcmp(Rtype{rr},'S11dB')
                result = invoke(mws,'Result1D','d1(1)1(1)');    % S11 in dB
                % Get nr of frequency points in the plot
                nRead = invoke(result,'GetN');
                [fin,S11in] = deal(zeros(nRead,1));
                for nn = 1:nRead
                    fin(nn) = invoke(result,'GetX',nn-1);        % Typically in GHz
                    S11in(nn) = invoke(result,'GetY',nn-1);
                end
                if isfield(M,'freq')
                    Nm = length(M.freq);
                    Rf{rr}.r = reshape(interp1(fin,S11in,M.freq,'spline'),Nm,1);
                    Rf{rr}.f = M.freq;
                else
                    Nm = nRead;
                    Rf{rr}.r = S11in;
                    Rf{rr}.f = fin;
                end
                Rf{rr}.t = Rtype{rr};
                release(result);
            elseif strcmp(Rtype{rr},'S11complex')
                resultA = invoke(mws,'Result1D','a1(1)1(1)');    % amplitude of S11
                resultP = invoke(mws,'Result1D','p1(1)1(1)');    % amplitude of S11
                % Get nr of frequency points in the plots
                nRead = invoke(resultA,'GetN');
                [fin,S11in] = deal(zeros(nRead,1));
                for nn = 1:nRead
                    fin(nn) = invoke(resultA,'GetX',nn-1);        % Typically in GHz
                    amp = invoke(resultA,'GetY',nn-1);
                    phase = rad(invoke(resultP,'GetY',nn-1));
                    S11in(nn) = amp.*exp(1i*phase);
                end
                if isfield(M,'freq')
                    Nm = length(M.freq);
                    Rreal = reshape(interp1(fin,real(S11in),M.freq,'spline'),Nm,1);
                    Rimag = reshape(interp1(fin,imag(S11in),M.freq,'spline'),Nm,1);
                    Rf{rr}.r = Rreal + 1i*Rimag;
                    Rf{rr}.f = M.freq;
                else
                    Nm = nRead;
                    Rf{rr}.r = S11in;
                    Rf{rr}.f = fin;
                end
                Rf{rr}.t = Rtype{rr};
                release(resultA);
                release(resultP);
            end
        end
        invoke(mws,'Save');
        invoke(mws,'Quit');
        
    case 'FEKO'
        % Build parameter string
        parStr = [];
        for nn = 1:Nn
            parStr = [parStr,' -#',M.params{nn},'=',num2str(xi(nn))];
        end
        % Remesh the structure with the new parameters
        FEKOmesh = ['cadfeko_batch ',[M.path,M.name,'.cfx'],parStr];
        system(FEKOmesh)
        % Run FEKO - cannot run with path, so change the directory
        curDir = pwd;
        cd(M.path)
        FEKOrun = ['runfeko ', [M.name,'.cfx']];
        system(FEKOrun)
        cd(curDir)
        % Generate output
        for rr = 1:Nr
            if strncmp(Rtype{rr},'S11',3)
                % Read the S11 touchstone file - must be exported by the FEKO
                % file with the correct name - Name_S11.s1p!
                [Spar,freq] = touchread([M.path,M.name,'_S11.s1p'],1);
                S11 = reshape(Spar(1,1,:),length(freq),1);
                Rf{rr}.f = freq;
            end
            if strcmp(Rtype{rr},'S11dB')
                Rf{rr}.r = dB20(S11);
            elseif strcmp(Rtype{rr},'S11complex')
                Rf{rr}.r = S11;
            end
            Rf{rr}.t = Rtype{rr};
        end
        
    case 'MATLAB'
        Ni = length(M.params);  % This is interpreted as the number of inputs to the function
        inType = [];
        for ii = 1:Ni
            inType = [inType,M.params{ii}];   % Initialise
        end
        switch inType
            case 'xf'
                Rfi = M.name(xi,M.freq);
            otherwise
                Rfi = M.name(xi);
        end
        % Distribute the responses
        % For MATLAB case the model should the return the specified responses
        % columnwise...
        for rr = 1:Nr
            Rf{rr}.r = Rfi(:,rr);
            Rf{rr}.t = Rtype{rr};
            if exist('f','var')
                Rf{rr}.f = f;
            end
        end
        
    otherwise
        error(['M.solver unknown for fine model evaluation'])
end

end

function Rc = coarseMod(M,xi,xp,f)

% Rc is a cell array of structures containing the response in Rc.r, the type Rc.t, and the
% (optional) domain (typically frequency) in Rc.f.  Same length as M.Rtype
% xi is an array of input parameters - same order as those specified in M
% xp is an array of implicit parameters - same order as those specified in M
% f is an array of frequency points where to evaluate the model (optional)
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'FEKO'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
%               'S11dB'
%               'S11complex'
%               'Gen'

% Limit the inputs - this should really never happen...
if isfield(M,'ximin')
    minI = xi < M.ximin;
    xi(minI) = M.ximin(minI);
    if minI | 0
		warning( strcat('Out of bounds coarse model evaluation encountered on ximin = ', ...
			mat2str(M.ximin), ', xi = ', mat2str(xi)) )
	end
end
if isfield(M,'ximax')
    maxI = xi > M.ximax;
    xi(maxI) = M.ximax(maxI);
    if maxI | 0
		warning( strcat('Out of bounds coarse model evaluation encountered on ximax', ...
			mat2str(M.ximax), ', xi = ', mat2str(xi)) )
	end
end
if isfield(M,'xpmin')
    minIp = xp < M.xpmin;
    xp(minIp) = M.xpmin(minIp);
    if minIp | 0
		warning( strcat('Out of bounds coarse model evaluation encountered on xpmin', ...
			mat2str(M.xpmin), ', xi = ', mat2str(xi)) )
	end
end
if isfield(M,'xpmax')
    maxIp = xp > M.xpmax;
    xp(maxIp) = M.xpmax(maxIp);
    if maxIp | 0
		warning( strcat('Out of bounds coarse model evaluation encountered on xpmax', ...
			mat2str(M.xpmax), ', xi = ', mat2str(xi)) )
	end
end

Nn = length(xi);
Nq = length(xp);
% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end
Nr = length(Rtype);
Rc = cell(1,Nr);

% Call the correct solver
switch M.solver
    case 'CST'
        error('CST solver not implemented yet for coarse model evaluations')
        
    case 'FEKO'
        % Build parameter string
        parStr = [];
        for nn = 1:Nn
            parStr = [parStr,' -#',M.params{nn},'=',num2str(xi(nn))];
        end
        % Also include possible implicit parameters
        for qq = 1:Nq
            parStr = [parStr,' -#',M.Iparams{qq},'=',num2str(xp(qq))];
        end
        
        % Remesh the structure with the new parameters
        FEKOmesh = ['cadfeko_batch ',[M.path,M.name,'.cfx'],parStr];
        system(FEKOmesh)
        % Run FEKO - cannot run with path, so change the directory
        curDir = pwd;
        cd(M.path)
        FEKOrun = ['runfeko ', [M.name,'.cfx']];
        system(FEKOrun)
        cd(curDir)
        % Generate output
        for rr = 1:Nr
            if strncmp(Rtype{rr},'S11',3)
                % Read the S11 touchstone file - must be exported by the FEKO
                % file with the correct name - Name_S11.s1p!
                [Spar,freq] = touchread([M.path,M.name,'_S11.s1p'],1);
                S11 = reshape(Spar(1,1,:),length(freq),1);
                Rc{rr}.f = freq;
            end
            if strcmp(Rtype{rr},'S11dB')
                Rc{rr}.r = dB20(S11);
            elseif strcmp(Rtype{rr},'S11complex')
                Rc{rr}.r = S11;
            end
            Rc{rr}.t = Rtype{rr};
        end
        
    case 'AWR'
        error('AWR solver not implemented yet for coarse model evaluations')
        
    case 'ADS'
        error('ADS solver not implemented yet for coarse model evaluations')
        
    case 'MATLAB'
        Ni = length(M.params);  % This is interpreted as the number of inputs to the function
        inType = [];
        for ii = 1:Ni
            inType = [inType,M.params{ii}];   % Initialise
        end
        switch inType
            case 'xxpf'
                Rci = M.name(xi,xp,f);
            case 'xf'
                Rci = M.name(xi,f);
            case 'xxp'
                Rci = M.name(xi,xp);
            otherwise
                Rci = M.name(xi);
        end
        % Distribute the responses
        % For MATLAB case the model should return the specified responses
        % columnwise...
        for rr = 1:Nr
            Rc{rr}.r = Rci(:,rr);
            Rc{rr}.t = Rtype{rr};
            if exist('f','var')
                Rc{rr}.f = f;
            end
        end
    otherwise
        error(['M.solver unknown for coarse model evaluation'])
end
end

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

end


