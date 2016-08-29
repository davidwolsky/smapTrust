function [Ri,Si,Pi,Ci,Oi,Li] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts)

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
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   The following (2) limits are only for warning generation - not used in optimization
%   ximin:  Vector of minimum values to limit the parameters [Nn,1] 
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% Mc:       Coarse model structure (can be cell array of structures if more than one type has to be calculated to get all the fine model responses)
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'FEKO'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
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
%   goalCent:   Cell array of centre point of goal domain {1,Ng} (used by the 'bw' goalType) (optional)
%   errNorm:    Cell array of type of error norms to use for optimization {1,Ng}
%               Valid types: (1,2,inf)
%   TolX:       Termination tolerance on X [ positive scalar - default 10^-2]
%   optsFminS:  options for the fminsearch local optimizer
%   optsPBIL:   options for the PBIL global optimizer
%   plotIter:   Flag to plot the responses after each iteration


% Returns:
% Ri:   Structure containing the responses at each iteration
% Si:   Structure containing the surrogates at each iteration
% Pi:   Structure containing the parameters at each iteration
% Ci:   Structure containing the costs at each iteration
% Oi:   Structure containing the 'optimizer' information at each iteration
% Li:   Structure containing the limiting information at each iteration

% Date created: 2015-03-06
% Dirk de Villiers and Ryno Beyers
% Last Modified: 2016-08-10
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
% 2016-08-21: Refactored the main loop to get fine model evaluation at the end and first iteration setup before loop.
%
% ToDo: TR
%       

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

if isfield(OPTopts,'Ni'), Ni = OPTopts.Ni; end
% TODO: DWW: 
if isfield(OPTopts,'TRNi'), TRNi = OPTopts.TRNi; end
if isfield(OPTopts,'TolX'), TolX = abs(OPTopts.TolX); end % Force positive
if isfield(OPTopts,'globOpt'), globOpt = OPTopts.globOpt; end
if isfield(OPTopts,'M_PBIL'), M_PBIL = OPTopts.M_PBIL; end
if isfield(OPTopts,'globOptSM'), globOptSM = OPTopts.globOptSM; end
if isfield(OPTopts,'optsFminS'), optsFminS = OPTopts.optsFminS; end
if isfield(OPTopts,'M_PBIL'), M_PBIL = OPTopts.M_PBIL; end
if isfield(OPTopts,'optsPBIL'), optsPBIL = OPTopts.optsPBIL; end
if isfield(OPTopts,'plotIter'), plotIter = OPTopts.plotIter; end


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
% 0) Normalise parameters
% 0) Optimise coarse model to find initial alignment position
% 0) Evaluate the fine model at starting position
% 0) Get the initial response
% 1) 4) Optimize the current model Si to find xi
% 2) 5) Test for convergence
% 3) 1) Evaluate the fine model at current position (Rfi)
% 4) 2) Get the response of the previous iteration surrogate at current
%    position (Rsi).  Only the position x(i-1) was calculated in the previous
%    iteration.
% 5) 3) Align the model at the current position the get Rsai and Si

%%% 1) Evaluate the fine model at current position (Rfi)
%%% 2) Get the response of the previous iteration surrogate at current
%%%    position (Rsi).  Only the position x(i-1) was calculated in the previous
%%%    iteration.
%%% 3) Align the model at the current position the get Rsai and Si
%%% 4) Optimize the current model Si to find xi
%%% 5) Test for convergence

% [Rci,Rfi,Rsi,Rsai,Si] = deal(cell(1,Ni));
% Rsai is the aligned surrogate, and Rsi the optimized surrogate at each
% iteration
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
Delta{1} = 0.25;
eta1 = 0.05;
eta2 = 0.9;
alp1 = 2.5;
alp2 = 0.25;

% Optimize coarse model to find initial alignment position
if globOpt
    [xinitn,costSi,exitFlag,output] = PBILreal(@(xin) costSurr(xin,Sinit,OPTopts),ximinn,ximaxn,M_PBIL,optsPBIL);
    xinitn = reshape(xinitn,Nn,1);
end
LHSmat = [];
RHSvect = [];
nonLcon = [];
[xin{1}, costS{1}] = fminsearchcon(@(xin) costSurr(xin,Sinit,OPTopts),xinitn,ximinn,ximaxn,LHSmat,RHSvect,nonLcon,optsFminS);
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

% Plot the initial fine, coarse, optimised surrogate and aligned surrogate
plotModels(plotIter, 1, Rci, Rfi, Rsi, Rsai, OPTopts);

% Test fine model response
costF{1} = costFunc(Rfi{1},OPTopts);
    
while ii <= Ni && ~specF && ~TolX_achieved
%Coming into this iteration as ii now with the fine model run here already and responses available. 

    % Exit if spec is reached (will typically not work for eq and never for minimax, and bw is explicitly excluded)
    % if costFi == 0 && isempty(find(ismember(OPTopts.goalType,'bw'),1))   
    if costF{ii} == 0 && isempty(find(ismember(OPTopts.goalType,'bw'),1))
        specF = 1;
    else
        specF = 0;
        
        % TR
        TRsuccess = 0;
        kk = 1;
        while ~TRsuccess && kk < TRNi
            % Set up TR boundaries
            ximinnTR = max((xin{ii} - Delta{ii}),ximinn);
            ximaxnTR = min((xin{ii} + Delta{ii}),ximaxn);
            
            % Do optimization
            % CRC_DDV: is this meant to always happen?
            if globOpt
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
            %         if TolX_achieved, keyboard; end
            
            enforceFineModelLimits();
            
            Rci{ii+1} = coarseMod(Mc,xi{ii+1},Sinit.xp,fc);
            Rfi{ii+1} = fineMod(Mf,xi{ii+1});
            
            for rr = 1:Nr
                % Get the surrogate response after previous iteration
                % optimization - thus at current iteration position
                Rsi{ii+1}{rr}.r = evalSurr(xi{ii+1},Si{ii+1-1}{rr});
                Rsi{ii+1}{rr}.t = Rci{ii+1}{rr}.t;
                if isfield(Rci{ii+1}{rr},'f'), Rsi{ii+1}{rr}.f = Rci{ii+1}{rr}.f; end
                if globOptSM < 2, SMopts.globOpt = 0; end
                if ~useAllFine
                    Si{ii+1}{rr} = buildSurr(xi{ii+1},Rfi{ii+1}{rr}.r,Si{ii+1-1}{rr},SMopts);
                else
                    for iii = 1:ii+1
                        r{iii} = Rfi{iii}{rr}.r;
                    end
                    Si{ii+1}{rr} = buildSurr(xi,r,Si{ii+1-1}{rr},SMopts);
                end
                % Also get the currently aligned surrogate for comparison
                Rsai{ii+1}{rr}.r = evalSurr(xi{ii+1},Si{ii+1}{rr});
                Rsai{ii+1}{rr}.t = Rci{ii+1}{rr}.t;
                if isfield(Rci{ii+1}{rr},'f'), Rsai{ii+1}{rr}.f = Rci{ii+1}{rr}.f; end
            end
            
            % Test fine model response
            costF{ii+1} = costFunc(Rfi{ii+1},OPTopts)
            costS{ii+1} = costSi
            
			% TODO: DWW: These values are far to high.
			%			 Also the first value is always wrongish i think
			%			 with cost of surrogate...
			
            % Evaluate results and adjust radius for next iteration
			% CRC_DDV: what about when both of these are negative?
			costChangeF = (costF{ii} - costF{ii+1})
			costChangeS = (costS{ii} - costS{ii+1})
			if ( costChangeF > 0 && costChangeS > 0 && abs(costChangeS) > TolX )
				rho{ii}{kk} = (costChangeF)./(costChangeS);
			else
				rho{ii}{kk} = 0.0;
			end
			
			rho{ii}
			Delta{ii}
            keyboard
			% TODO: DWW: These values seem too low.
            sk{ii} = xin{ii+1}-xin{ii};
            if rho{ii}{kk} >= eta2
                TRsuccess = 1;
                Delta{ii+1} = max(alp1.*norm(sk{ii}),Delta{ii})
            elseif rho{ii}{kk} > eta1
                TRsuccess = 1;
                Delta{ii+1} = Delta{ii}
            else
                TRsuccess = 0;
                Delta{ii} = alp2.*norm(sk{ii}) % Shrink current Delta
            end
            kk = kk+1;
        end
        
        % Make a (crude) log file
        %     save SMlog ii xi Rci Rfi Rsi Si costS costF limF limC limS
        save SMlog ii xi Rci Rfi Rsi Si costS costF limMin_f limMax_f limMin_c limMax_c rho  Delta
        
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

Ci.costS = costS;
Ci.costF = costF;

Oi.specF = specF;
Oi.TolX_achieved = TolX_achieved;   % Flag
Oi.TolXnorm = TolXnorm; % Actual value
Oi.Ni = ii;
Oi.rho = rho;
Oi.Delta = Delta;


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


function plotModels(plotFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts)
    if plotFlag
        figure(itNum)
		Nr = length(OPTopts.Rtype); % Number of responses requested
        for rr = 1:Nr
            subplot(Nr,1,rr)
            if isfield(Rci{itNum}{rr},'f')
                plot(Rci{itNum}{rr}.f,Rfi{itNum}{rr}.r,'k'), grid on, hold on
                plot(Rci{itNum}{rr}.f,Rci{itNum}{rr}.r,'r'), grid on, hold on
                plot(Rsi{itNum}{rr}.f,Rsi{itNum}{rr}.r,'b')
                plot(Rsai{itNum}{rr}.f,Rsai{itNum}{rr}.r,'g--')
                xlabel('Frequency')
                % Plot the specs...
            else
                plot(Rfi{itNum}{rr}.r,'k'), grid on, hold on
                plot(Rci{itNum}{rr}.r,'r'), grid on, hold on
                plot(Rsi{itNum}{rr}.r,'b')
                plot(Rsai{itNum}{rr}.r,'g')
                xlabel('Index')
                % Plot the specs...
            end
            ylabel(OPTopts.Rtype{rr})
            title(['Iteration ',num2str(itNum)])
            legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate')
        end

    end
end % plotModels


function Rf = fineMod(M,xi)

% Rf is a cell array of structures containing the response in Rf.r, the type Rf.t, and the
% (optional) domain (typically frequency) in Rf.f.  Same length as M.Rtype
% xi is an array of input parameters - same order as those specified in M
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB'/'FEKO' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
%               'S11dB' - obvious!

% Limit the inputs
if isfield(M,'ximin')
    minI = xi < M.ximin;
    xi(minI) = M.ximin(minI);
    if minI | 0
		warning( strcat('Out of bounds fine model evaluation encountered on ximin = ', ...
			mat2str(M.ximin), ', xi = ', mat2str(xi)) )
	end
	
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


