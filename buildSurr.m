function Si = buildSurr(xi,Rfi,S,opts)

% function Si = buildSurr(xi,Rfi,S,opts)
% function to build a surrogate model structure S, using the input fine 
% model response pairs {xi,Rfi}.  More than one input response pair can be
% specified, as is typical in space mapping modelling instead of
% optimization applications.  User specified options - mainly flags to
% specify SM types - are provided in opts.
%
% Based on the 2006 MTT paper by Koziel et al on a SM Framework
%
% In general the extraction performs best when all input variables,
% including implicit variables, are of the degree 1.  This is because the
% internal optimizer attempts to optimize additive as well as multiplicative
% variables - the latter typically normalized to 1 - at the same time.
%
% Inputs:
% xi: Response positions {1,Ni}[Nn,1] (Nn input parameters, Ni previous fine evaluations)
% Rfi: Fine model responses {1,Ni}[Nm,1] (Nm output responses, Ni previous fine evaluations)
%   xi and Rfi may be cell arrays of the same length indicating multiple fine
%   model evaluations.  The surrogate will then be a best fit over the entire
%   extended parameter space.
% S:  Initial surrogate model containing at least
%  S.coarse -  Function handle pointing to a predefined function that evaluates
%          the coarse model operating as Rc = coarse(xc,xp,f) thus returning Rc 
%          (coarse model response) for each element in xc. Optional xp are
%          the implicit (pre-assigned) parameters.
%  S.xp - Optional pre-assigned variables (depending on the coarse model [Nq,1]
%  S.f - Optional frequency vector where the responses are calculated [Nm,1] (must be included for FSM)
% opts: User options structure - mostly flags to specify type of SM (defaults)
%   getA - flag to get multiplicative OSM (set to 2 to get M factors), typically 1 to force 
%          single factor for all outputs which is much faster and often works well (0).
%   getB - flag to get multiplicative input SM (0)
%          If getB == 1 the full matrix will be extracted
%          If getB == 2 only diagonal entries will be extracted (uncoupled system)
%          If getB is a boolean vector of length Nn, only the 'true'
%          entries on the diagonal will be extracted.  This is an uncoupled
%          system with only a subset of the parameters allowed to have
%          a linear variation.
%   getc - flag to get additive input SM (1)
%   getG - flag to get multiplicative part of ISM (0)
%          If getG == 1 the full matrix will be extracted
%          If getG is a boolean vector of length Nn, only the 'true'
%          columns be extracted.  This can be used to remove dependence on
%          certain of the input parameters. 
%          If getG is a boolean Matrix of size [Nq,Nn], only the 'true'
%          entries will be optimized.  This allow the user full control of
%          dependencies between optimization and implicit parameters
%   getxp- flag to optimize additive part of ISM (0) 
%   getd - flag to get additive OSM (0)
%   getE - flag to get first order OSM (0)
%   getF - flag to get frequency mapping (0)
%   wk - weights to determine the error function in the model fitting when
%        more than one fine point is included. Can be the same length 
%        as the number of cells in xi and Rfi. Default all (1). See eq
%        (11) in the reference. If length(xi) > 1 and length(wk) == 1 then
%        wk will be generated internally as an exponential function,
%        growing by the factor wk ie: wk.^[1:length(xi)]
%        Additional option added for wk == 0 with length(xi) == 1. This option
%        uses the number of SM types requested (implying the number of unknowns)
%        to keep track of how many fine model points to include. Ones are used for
%        the weight of the number of SM types for the most recent fine points 
%        and a zero weight is applied to all previous entries.
%   vk - weights to determine the error function in the Jacobian fitting when
%        more than one fine point is included. Can be the same length 
%        as the number of cells in xi and Rfi. Default all (0). See eq
%        (11) in the reference. If length(xi) > 1 and length(vk) == 1 then
%        vk will be generated internally as an exponential function,
%        growing by the factor vk ie: vk.^[1:length(xi)]
%   ximin - vector of minimum values for xi to constrain the search space [Nn,1]
%   ximax - vector of maximum values for xi to constrain the search space [Nn,1]
%   xpmin - vector of minimum values for xp to constrain the search space [Nn,1]
%   xpmax - vector of maximum values for xp to constrain the search space [Nn,1]
%           All the above limits are on the coarse model space - they are
%           shifted internally for the surrogate model space.  That is to
%           say the coarse model will only ever be evaluated within these
%           bounds during the PE process.  Equal upper and lower bounds
%           will remove the parameters of interest from the PE optimization
%   Amin - minimum value for A to constrain the search space (default 0.5)
%   Amax - maximum value for A to constrain the search space (default 2.0)
%          The above values can be scalar of vectors of length Nm
%   Fmin - minimum values for F to constrain the search space [2,1]
%   Fmax - maximum values for F to constrain the search space [2,1]
%          Once again, equal upper and lower bounds can be used above to
%          remove a certain element from the search space
%   fmin - minimum shifted frequency when using getF
%   fmax - maximum shifted frequency when using getF
%   optsFminS - options for the fminsearch local optimizer in the PE
%   globOpt - flag to include global search (PBILreal) in the PE (Not currently working)
%   M_PBIL - number of bits per parameter in the PBIL global search (8)
%   optsPBIL - options for the PBIL global optimizer in the PE
%   errNorm - Type of error norm to use for parameter extraction (1,2,inf)
%   errW - Vector of weights (typically binary but can be any real number),
%          of length Nm, to calculate the extraction error.  Can be used to 
%          mask out regions in the response domain of less importance. 
%          Default ones. 
%   plotAlignmentFlag - When true a plot of the starting point and final point
%           of alignment will be plotted. For complex results both the real and 
%           imaginary part is shown along with the error between the fine model and
%           the evaluated surrogate. NOTE an extra two coarse model evaluations are 
%           done when populating the errors to plot.
%
% Returns:
% S -  The surrogate model structure containing:
% coarse:  Function handle pointing to a predefined function that evaluates
%          the coarse model operating as Rc = coarse(xc,xp) thus returning Rc 
%          (coarse model response) for each element in xc. Optional xp are
%          the implicit (pre-assigned) parameters.
% A:       Multiplicative OSM factor diag[Nm,Nm]
% B:       Multiplicative input SM factor [Nn,Nn]
% c:       Additive input SM term [Nn,1]
% G:       Multiplicative ISM factor [Nq,Nm]
% xp:      Pre-assigned parameters [Nq,1]
% d:       Additive zeroth order OSM term [Nm,1]
% E:       Additive first order OSM term [Nm,Nn]
% xi:      Position of last update point of the model [Nn,1]
% F:       Frequency space mapping parameters [2,1]
% freq:    Coarse model frequency range [Nm,1]
%
% Date created: 2014-11-09
% Dirk de Villiers
% Last Modified: 2016-03-07
% Updates:
% 2014-11-09: Write function shell and basic functionality
% 2014-11-10: Fix constraints of input SM (case 1), cases 2 - 5
% 2014-11-11: Implement add OSM (d), Fix B to diagonal matrix
% 2014-11-13: Add some more 2 variable cases. Include global search flag
%             and constraint flags. Add global search to case 5.
% 2014-11-14: Add error norm and error weight options
% 2014-11-22: Move some general code out of switch
% 2014-11-23: Include optimizer options as inputs
% 2015-02-23: Fix implicit constraints
% 2015-03-03: Fix simple special case of only getA for one fine model
%             Include c and B cases, and simplify Bc case (no constraints)   
% 2015-03-04: Include AB, Ac, ABc, G, xp, cxp cases
%             Simplify Gxp case (no constraints)
%             Include AlimMin/Max as opts parameter 
%             Include Gv_init variable for consistency with Bv_init
% 2015-03-05: Include globalOpt in BcGxp case
%             Include BG, cG, BcG, BGxp, cGxp cases
%             REDESIGN CASES COMPLETELY.  ALL GENERAL.
% 2015-03-09: Start FSM
% 2015-03-10: Continue FSM
% 2015-03-11: Finish FSM
% 2015-03-12: Constrained search added
% 2015-03-13: Single value wk functionality added
% 2015-03-14: Fix errW NaN issue, remove limits from parameters space-
%             should be handled externally
% 2015-04-04: Fix matrix size issue when getA == 1 and more than one
%             iteration is required. A then becomes a scalar, and Nm was
%             calculated as 1. Have to run coarse model every time.
% 2015-05-07: Limit xi and xp in the error function (based on Alex
%             Vermeulen's code)
% 2015-10-16: Fixed small issue with xplims when none are provided
% 2016-02-24: Fixed small issue with xplims when no xp is provided
%             Allow off-diagonal B entries    
% 2016-02-26: Fixed serious issue with FSM
% 2016-02-29: Start redesign to include the bounds while using
%             fminsearchcon
% 2016-03-01: Continue redesign to handle bounds properly - get rid of
%             globOpt option for now
% 2016-03-07: Continue with new bounds design - change standard workflow to
%             always include all the types of SM in the S model.  Keep defaults if not
%             requested in PE.
% 2017-03-08: Introduced a wk == 0  option that used the number of SM types to determine how
%             how many of the most recent fine points to use when calculating the error function
%             for model fitting. 
% 2017-07-08: Bring in functions from the Matlab Optimisation toolbox. To do this the constraints
%             and variables need to take a different form. 
% 2017-07-29: Removed fixed parameters from the optimisation space. This is especially neccessary for
%             the global optimiser and to reduce the complexity for debugging.
% 2017-08-14: Introduced a plotting flag and function to make a graph showing the error before and 
%             after the alignment phase.

% ToDo: Implement E (first order OSM)
% ToDo: Jacobian fitting in error functions (vk)
% ToDo: Re-introduce globOpt option
% ToDo: Include more optimizer options - start with fmincon

% Preassign some variables
[lenA,lenB,lenc,lenG,lenxp,lenF] = deal(0);
[Nc,Nn,Nq,Nm] = deal(0);
[A_init,Av_init,B_init,Bv_init,c_init,G_init,Gv_init,xp_init,F_init,typeA,typeB,typec,typeG,typexp,typeF] = deal([]);
[Bmin,Bmax,cmin,cmax,Gmin,Gmax,xpmin,xpmax,fmin,fmax] = deal([]);

% Sort out formats
if iscell(xi) && iscell(Rfi) && (length(xi) == length(Rfi))   % Basic error checking - if any issue here only first point will be used
    Nc = length(xi);  % Number of input point cells
else     % Force only first point to be used, and make cell arrays
    if ~iscell(xi), xi = {xi}; end 
    if ~iscell(Rfi), Rfi = {Rfi}; end 
    Nc = 1;
end

% Get vector sizes
Nn = length(xi{Nc});    % Number of input parameters
Nm = length(Rfi{1});    % Number of responses

% Default SM type requests
getA = 0;
getB = 0;
getc = 1;   
getG = 0;
getxp = 0;
getd = 0;
getE = 0;
getF = 0;

plotOpts = {}

% Default constraints
ximin = -inf.*ones(Nn,1); 
ximax = inf.*ones(Nn,1);
if isfield(S,'xp')
    if ~isempty(S.xp)
        Nq = length(S.xp);
        xpmin = -inf.*ones(Nq,1);
        xpmax = inf.*ones(Nq,1);
    else isempty(S.xp)
        S = rmfield(S,'xp');    % Get rid of empty xp field for later checks
    end
end
Amin = 0.5;
Amax = 2.0;

% Get user SM type request
if isfield(opts,'getA'), getA = opts.getA; end
if isfield(opts,'getB'), getB = opts.getB; end
if isfield(opts,'getc'), getc = opts.getc; end
if isfield(opts,'getG'), getG = opts.getG; end
if isfield(opts,'getxp'), getxp = opts.getxp; end
if isfield(opts,'getd'), getd = opts.getd; end
if isfield(opts,'getE'), getE = opts.getE; end
if isfield(opts,'getF'), getF = opts.getF; end

if isfield(S,'f')
    % TODO_DWW: Move to somewhere more logical or bac up and fix usage logic with getF
    if getF
        fmin = 0.9*min(S.f);
        fmax = 1.1*max(S.f);
    end
end

NSMUnknowns = getNSMUnknowns();

% Default optimization parameters
optsFminS = optimset('display','none');
globOpt = 0;
% M_PBIL = 8;
% optsPBIL = [];
% CRC_DDV: I dont think 2 is the correct default?
errNorm = 2;
errW = 1;
plotAlignmentFlag = 0;
if isfield(opts,'wk') 
    wk = opts.wk
    if length(wk) == 1
        if (wk == 0) 
            wk = ones(1,Nc);
            % If the error function for the model fitting becomes overdetermined 
            % and can give a skewed model.
            if (Nc > NSMUnknowns)
                wk(1:end-NSMUnknowns) = 0
            end
        else
            wk = wk.^[1:Nc]
        end
    end
else
    wk = ones(1,Nc)
end

if isfield(opts,'vk') 
    % CRC_DDV: DWW: Should the weighting for model fitting be applied to this Jacobian fitting too?
    vk = opts.vk;
    if length(vk) == 1
        vk = vk.^[1:Nc];
    end
else
    % CRC_DDV: DWW: The documentation says that this should be zeros by default?
    vk = ones(1,Nc);
end

% Get user constraints
if isfield(opts,'ximin'), ximin = opts.ximin; end
if isfield(opts,'ximax'), ximax = opts.ximax; end
if isfield(opts,'Amin'), Amin = opts.Amin; end
if isfield(opts,'Amax'), Amax = opts.Amax; end
if isfield(opts,'xpmin'), xpmin = opts.xpmin; end
if isfield(opts,'xpmax'), xpmax = opts.xpmax; end
if isfield(opts,'Fmin'), Fmin = opts.Fmin; end
if isfield(opts,'Fmax'), Fmax = opts.Fmax; end
% TODO_DWW: Move this to inside an getF scope
% if isfield(opts,'fmin'), fmin = opts.fmin; end
% if isfield(opts,'fmax'), fmax = opts.fmax; end

% Get user optimization parameters
if isfield(opts,'optsFminS'), optsFminS = opts.optsFminS; end
if isfield(opts,'globOpt'), globOpt = opts.globOpt; end
% if isfield(opts,'M_PBIL'), M_PBIL = opts.M_PBIL; end
% if isfield(opts,'optsPBIL'), optsPBIL = opts.optsPBIL; end
if isfield(opts,'errNorm'), errNorm = opts.errNorm; end
if isfield(opts,'errW'), errW = opts.errW; end
if isfield(opts,'plotAlignmentFlag'), plotAlignmentFlag = opts.plotAlignmentFlag; end

% Assign current iteration S-model
Si = S;

% Initialise the optimization variable vector
% The limits named *min/*max indicate the global search range.  Actual
% values may fall outside these ranges after the local search alignment...
if getA
    typeA = 'A';
    if getA == 1 % Use a single A for all responses
        lenA = 1;
        A_init = 1;
        if isfield(S,'A')
            A_init = S.A;
        end
        Av_init = mean(diag(A_init));   % In case a matrix was provided
        A_init = Av_init;               % In case a matrix was provided
        Amin = min(Amin);
        Amax = max(Amax);
    elseif getA == 2
%         % Have to run the coarse model to find the response size
%         Rc = evalSurr(xi{Nc},S);
%         [Nm,Np] = size(Rc);
%   want to take the number of paramaters requried as the number of actual values used. 
        lenA = Nm;
        if isfield(S,'A')
            A_init = S.A;
            Av_init = diag(A_init);
            if length(Av_init) ~= 1 && length(Av_init) ~= Nm
                error('Size of provided S.A not compatible with response size');
            elseif length(Av_init) == 1     % Expand to right length
                Av_init = ones(Nm,1).*Av_init;
                A_init = diag(Av_init);
            end
        else
            A_init = eye(Nm);
            Av_init = diag(A_init);
        end
        % Expand the limits if needed - user supplied A overrides all else
        if length(Amin) ~= 1 && length(Amin) ~= Nm
            warning('Length of supplied Amin incompatible with length of requested A - minimum value used');
            Amin = ones(Nm,1).*min(Amin);
        end
        if length(Amax) ~= 1 && length(Amax) ~= Nm
            warning('Length of supplied Amax incompatible with length of requested A - maximum value used');
            Amax = ones(Nm,1).*min(Amax);
        end
        % Expand the limits if they are scalar
        if length(Amin) == 1
            Amin = ones(Nm,1).*Amin;
        end
        if length(Amax) == 1
            Amin = ones(Nm,1).*Amax;
        end
    else
        error(['Unknown getA flag: ', num2str(getA),', should be 0, 1 or 2']);
    end
else    % A not requested - revert to default and clamp box limits.  Do here since user can supply limits
    [A_init,Av_init,Amin,Amax] = deal(1);
    lenA = 1;
end

% --- getB ----
lenB = Nn*Nn;
B_init = eye(Nn);
[Bv_init,Bmin,Bmax] = deal(reshape(B_init,lenB,1));
BminDefaultDiag = 0.5;
BmaxDefaultDiag = 2.0;
BminDefaultCross = -0.5;
BmaxDefaultCross = 0.5;
if any(getB)
    typeB = 'B';
    if isfield(S,'B')
        B_init = S.B;
    end
    % Sort out the box limits according to what is requested
    if numel(getB) == 1 && getB == 1    % Full matrix to be optimized
        BminM = BminDefaultCross.*ones(Nn,Nn).*~eye(Nn) + diag(BminDefaultDiag.*ones(Nn,1));  
        BmaxM = BmaxDefaultCross.*ones(Nn,Nn).*~eye(Nn) + diag(BmaxDefaultDiag.*ones(Nn,1));
    elseif numel(getB) == 1 && getB == 2    % Only the diagonal entries optimized
        % Keep the rest of the entries the same as B - this allows (rarely
        % used) the off-diagonal entries to be the user supplied values
        [BminM,BmaxM] = deal(B_init);
        BminM = BminM.*~eye(Nn) + diag(BminDefaultDiag.*ones(Nn,1));
        BmaxM = BmaxM.*~eye(Nn) + diag(BmaxDefaultDiag.*ones(Nn,1));
    elseif numel(getB) == Nn   % Only certain of the diagonal entries optimized
        % Keep the rest of the entries the same as B - this allows (rarely
        % used) the off-diagonal entries to be the user supplied values
        [BminM,BmaxM] = deal(B_init);
        [minDiag,maxDiag] = deal(zeros(Nn,1));
        minDiag(getB == 1) = BminDefaultDiag;
        maxDiag(getB == 1) = BmaxDefaultDiag;
        BminM = BminM.*(~diag(getB)) + diag(minDiag);
        BmaxM = BmaxM.*(~diag(getB)) + diag(maxDiag);
    else
        error(['Unknown getB flag: ', num2str(getB),', should be 0, 1, 2 or a bollean vector of length Nn = ', num2str(Nn)]);
    end

    % Override defaults with user input
    if isfield(opts,'Bmin'), BminM = opts.Bmin; end
    if isfield(opts,'Bmax'), BmaxM = opts.Bmax; end

    Bv_init = reshape(B_init,lenB,1);
    Bmin = reshape(BminM,lenB,1);
    Bmax = reshape(BmaxM,lenB,1);
end

% --- getC ----
lenc = Nn; 
[c_init,cmin,cmax] = deal(zeros(lenc,1));
if getc
    typec = 'c';
    if isfield(S,'c')
        c_init(:,1) = S.c;
    end

    cmin = ximin - reshape(Bmax, Nn, Nn)*ximax;
    cmax = ximax - reshape(Bmin, Nn, Nn)*ximin;

    % Override defaults with user input
    if isfield(opts,'cmin'), cmin = opts.cmin; end
    if isfield(opts,'cmax'), cmax = opts.cmax; end
end

% --- getxp and getG ----
% Have to provide an xp input for implicit space mapping...
if isfield(S,'xp')
    lenxp = Nq;
    % Note name pmin/max here - xpmin reserved for linear constraints
    [xp_init(:,1),pmin(:,1),pmax(:,1)] = deal(S.xp);    
    lenG = Nq*Nn;
    G_init = zeros(Nq,Nn); 
    [Gv_init,Gmin,Gmax] = deal(reshape(G_init,lenG,1));

    GminDefault = -2.0;
    GmaxDefault = 2.0;
    if any(getG)
        typeG = 'G';
        if isfield(S,'G')
            G_init = S.G;
        end
        % Sort out the box limits according to what is requested
        if numel(getG) == 1 && getG == 1    % Full matrix to be optimized
            GminM = GminDefault.*ones(Nq,Nn);
            GmaxM = GmaxDefault.*ones(Nq,Nn);
        elseif numel(getG) == Nn  % Only certain input parameter dependencies are kept for all implicit parameters
            % Keep the rest of the entries the same as G
            [GminM,GmaxM] = deal(G_init);
            % Unconstrained columns where we want to allow a search
            GminM(:,getG == 1) = GminDefault;
            GmaxM(:,getG == 1) = GmaxDefault;
        elseif all(size(getG) == [Nq,Nn])
            % Keep the rest of the entries the same as G
            [GminM,GmaxM] = deal(G_init);
            % Unconstrained entries where we want to allow a search
            GminM(getG == 1) = GminDefault;
            GmaxM(getG == 1) = GmaxDefault;
        else
            error(['Unknown getG flag: ', num2str(getG),', should be 0, 1 or a bollean vector of length Nn = ', num2str(Nn),', or a boolean matrix of size [Nq,Nn] = [',num2str(Nq),',',num2str(Nn),']']);
        end
        
        % Override defaults with user input
        if isfield(opts,'Gmin'), GminM = opts.Gmin; end
        if isfield(opts,'Gmax'), GmaxM = opts.Gmax; end
        
        Gv_init = reshape(G_init,lenG,1);
        Gmin = reshape(GminM,lenG,1);
        Gmax = reshape(GmaxM,lenG,1);
    end

    if getxp
        typexp = 'xp';
        xp_init(:,1) = S.xp;
        % Box constraint
        pmin = xpmin - reshape(Gmax, Nq, Nn)*ximax;
        pmax = xpmax - reshape(Gmin, Nq, Nn)*ximin;
    end

else
    if getxp
        warning('Cannot perform getxp with no S.xp provided - request ignored');
    elseif any(getG)
        warning('Cannot perform getG with no S.xp provided - request ignored');
    end
    [xp_init,pmin,pmax,G_init,Gv_init,Gmin,Gmax] = deal([]);
end

% --- getF ---
lenF = 2;
if getF 
    if ~isfield(S,'f')
        error('If getF == 1 you have to supply a frequency vector in the SM model');
    end
    typeF = 'F';
    if isfield(S,'F')
        F_init = S.F;
    else
        F_init = [1;0];
    end
    F1min = -0.5;
    F1max = 0.5;
    F2min = fmin - F1max*max(S.f);
    F2max = fmax - F1min*min(S.f);
    Fmin = [F1min; F2min];
    Fmax = [F1max; F2max];

else  % F not requested - revert to default and clamp box limits.  Do here since user can supply limits.
    [F_init,Fmin,Fmax] = deal([1;0]);
end

% Get positions of different parts of the optimization vector
lenVect = [lenA,lenB,lenc,lenG,lenxp,lenF];
[firstPos,lastPos] = deal(zeros(1,6));
for pp = 1:6 
%     if pp > 1
        firstPos(pp) = 1 + sum(lenVect(1:(pp-1)));
%     end
    lastPos(pp) = sum(lenVect(1:pp));
end

initVect = [Av_init;Bv_init;c_init;Gv_init;xp_init;F_init];
minVect = [Amin;Bmin;cmin;Gmin;pmin;Fmin];
maxVect = [Amax;Bmax;cmax;Gmax;pmax;Fmax];
inputType = [typeA,typeB,typec,typeG,typexp,typeF];
Si.inputType = inputType;

% Set up the general parameter extraction options
optsParE.Nn = Nn;
optsParE.Nq = Nq;
optsParE.lenVect = lenVect;
optsParE.firstPos = firstPos;
optsParE.lastPos = lastPos;
optsParE.errNorm = errNorm;
optsParE.errW = errW;

% Set up and run the optimizations (parameter extractions)
% Use provided SM parameters as initial values if available, and enhance
% with global search if required
if strcmp(inputType,'F') || strcmp(inputType,'AF') && Nc == 1 % Special cases where the coarse model is not re-evaluated (for each optimization iteration).  Use interpolation/extrapolation instead...
    % Calculate the coarse model - make sure F is set to the default!
    % If not it will be shifted twice...
    S.F = [1,0];
    Rc = evalSurr(xi{Nc},S);
    
    % if 0 && globOpt    % Only use local optimizer for FSM
        % F_init = PBILreal(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),Fmin,Fmax,M_PBIL,optsPBIL);
    % end
%     Fvect = fminsearchcon(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),F_init,[0, -inf],[],[-min(S.f),0;0,-1],[0;0],[],optsFminS);     % Positive multiplier, and minimum frequency
    % TODO_DWW: Follow through here...
    Fvect = fminsearchcon(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),F_init,[0, -inf],[],[-min(S.f),-1],[0],[],optsFminS);     % Positive multiplier, and minimum frequency
    
    fs = Fvect(1).*S.f + Fvect(2);
    RsComp = interp1(S.f,Rc,fs,'linear');
    RsComp = reshape(RsComp,length(RsComp),1);
    % Get rid of NaNs from shift
    cleanPos = find(~isnan(RsComp));
    RsCompClean = RsComp(cleanPos);
    fClean = S.f(cleanPos);
    fClean = reshape(fClean,length(fClean),1);
    % Include the edge points for the interpolation
    if max(fClean) < max(S.f)
        RsCompClean = [RsCompClean;Rfi{1}(end)];
        fClean = [fClean;S.f(end)];
    end
    if min(fClean) > min(S.f)
        RsCompClean = [Rfi{1}(1);RsCompClean];
        fClean = [S.f(1);fClean];
    end
    Rs = interp1(fClean,RsCompClean,S.f,'spline');

    optVect = initVect;
    if strcmp(inputType,'AF')       % Need to also get the A factor
        A = Rfi{1}./reshape(Rs,Nm,1);
        if getA == 1, A = mean(A); end
        optVect(firstPos(1):lastPos(1)) = A;
    end
    optVect(firstPos(6):lastPos(6)) = reshape(Fvect,lenF,1);
elseif strcmp(inputType,'A') && Nc == 1  % Special case without optimization
    if ~exist('Rc','var')    % Check if coarse model has been calculated
        Rc = evalSurr(xi{Nc},S);
    end
    A = Rfi{1}./Rc;
    if getA == 1, A = mean(A); end
    optVect = initVect;
    optVect(firstPos(1):lastPos(1)) = A;
else
    % Set up the linear constraints
    keyboard
    Ncon = 2*Nn + 2*Nq + 2*(getF>0);
    LHS_mat = zeros(Ncon,length(initVect));
    
    lhsA_mat = zeros(size(A_init));        
    lhsB_mat = zeros(size(B_init));    % Basic matrix shape to use for distribution of the input vector in the LHS matrix for B
    lhsc_mat = zeros(size(c_init));    % Basic matrix shape to use for distribution of the input vector in the LHS matrix for c
    lhsG_mat = zeros(size(G_init));    % Basic matrix shape to use for distribution of the input vector in the LHS matrix for G
    lhsxp_mat = zeros(size(xp_init));
    lhsF_mat = zeros(size(F_init));
    
    RHS_vect = [-ximin; ximax; -xpmin; xpmax; -fmin; fmax];
    
    % First populate the input space limits in LHS_mat
    for nn = 1:Nn
        xA_vect = diag(lhsA_mat);       % Always zeros - no influence on input/implicit space or frequency bounds
        
        xB_mat = lhsB_mat;
        if lenB > 0;
            xB_mat(nn,:) = xi{1}';
        end
        xB_vect = reshape(xB_mat,1,lenB);
        
        xc_mat = lhsc_mat;
        if lenc > 0
            xc_mat(nn,:) = 1;
        end
        xc_vect = reshape(xc_mat,1,lenc);
        
        xG_vect = reshape(lhsG_mat,1,lenG);
        xxp_vect = reshape(lhsxp_mat,1,lenxp);
        xF_vect = reshape(lhsF_mat,1,lenF);
        
        xLBrow = [xA_vect,-xB_vect,-xc_vect,xG_vect,xxp_vect,xF_vect];
        xUBrow = [xA_vect,xB_vect,xc_vect,xG_vect,xxp_vect,xF_vect];
        
        LHS_mat(nn,:) = xLBrow; % Lower bound row
        LHS_mat(Nn+nn,:) = xUBrow; % Upper bound row
    end
    % And now the implicit parameters part
    for qq = 1:Nq
        % Always zeros
        xpA_vect = diag(lhsA_mat);
        xpB_vect = reshape(lhsB_mat,1,lenB);
        xpc_vect = reshape(lhsc_mat,1,lenc);
        xpF_vect = reshape(lhsF_mat,1,lenF);
        
        xpG_mat = lhsG_mat;
        if lenG > 0
            xpG_mat(qq,:) = xi{1}';
        end
        xpG_vect = reshape(xpG_mat,1,lenG);
        
        xpxp_mat = lhsxp_mat;
        if lenxp > 0
            xpxp_mat(qq,:) = 1;
        end
        xpxp_vect = reshape(xpxp_mat,1,lenxp);
                
        xpLBrow = [xpA_vect,xpB_vect,xpc_vect,-xpG_vect,-xpxp_vect,xpF_vect];
        xpUBrow = [xpA_vect,xpB_vect,xpc_vect,xpG_vect,xpxp_vect,xpF_vect];

        rLB = 2*Nn + qq;    % Lower bound row
        rUB = rLB + Nq;     % Upper bound row
        
        LHS_mat(rLB,:) = xpLBrow; % Lower bound row
        LHS_mat(rUB,:) = xpUBrow; % Upper bound row
    end
    % And frequency shifts
    if getF
        LHS_mat(2*(Nn+Nq)+1,end-1:end) = [-min(S.f),-1];
        LHS_mat(2*(Nn+Nq)+2,end-1:end) = [max(S.f),1];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display(['--- TODO_DWW: Starting parameter extraction ---'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO_DWW: Normalisation options from FEKO chat:
%   - normalise internal to the optimisation algorithms with standard deviation.
%   - normalise incoming between 0 and 1.
%   - the goal functions are also normalised.
    
    keyboard
    fullProblem = {};
    fullProblem.x0 = initVect;
    fullProblem.Aineq = LHS_mat;
    fullProblem.bineq = RHS_vect;
    fullProblem.lb = minVect;
    fullProblem.ub = maxVect;

    % TODO_DWW: Clean up
    [nProblem] = normaliseProblem(fullProblem, optsParE);
%     [x0] = denormaliseProblem(nProblem, fullProblem, optsParE);

    [reducedProblem] = removeFixedParameters(nProblem);


    % TODO_DWW: Clean up - not used if not plotted... 
    startingError = 0;
    completionError = 0;
    if ( plotAlignmentFlag == 1 )
        % Plot initial error before alignment
        plotOpts.plotTitle = 'Starting alignment';
        startingError = erri(reducedProblem.x0,xi,Rfi,S,wk,vk,optsParE, fullProblem, plotAlignmentFlag, plotOpts)
    end

    problem = {};
    problem.x0 = reducedProblem.x0;
    problem.Aineq = reducedProblem.Aineq;
    problem.bineq = reducedProblem.bineq;
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = reducedProblem.lb;
    problem.ub = reducedProblem.ub;
    problem.nonlcon = [];
    if globOpt
        problem.solver = 'ga';
        problem.fitnessfcn = @(tempOptVect) erri(tempOptVect,xi,Rfi,S,wk,vk,optsParE,fullProblem, false, plotOpts);
        problem.nvars = length(reducedProblem.x0);
        problem.options = optimoptions('ga')
        [optVectGlobalReduced,fval,exitflag,output] = ga(problem)
        % Start with global search to get initial value
        problem.x0 = optVectGlobalReduced
    end
    problem.solver = 'fmincon';
    problem.objective = @(tempOptVect) erri(tempOptVect,xi,Rfi,S,wk,vk,optsParE, fullProblem, false, plotOpts);
    problem.options = optimoptions('fmincon',...
                                   'Display','iter-detailed',...
                                   'Diagnostics','on',...
                                   'PlotFcn',@optimplotconstrviolation)
    problem.options.DiffMinChange = 1e-4;
%     problem.options.Algorithm = 'sqp';
    % keyboard
    [optVectReduced,fval,exitflag,output] = fmincon(problem)
    if ( plotAlignmentFlag == 1 )
        % Plot errors after alignment
        plotOpts.plotTitle = 'Alignment complete'
        completionError = erri(optVectReduced,xi,Rfi,S,wk,vk,optsParE, fullProblem, plotAlignmentFlag, plotOpts);

        assert(completionError == fval, 'The optimised parameters should return the same error as the optimiser received.')
        % TODO_DWW: Clean up - not used if not plotted... 
%         assert(startingError > completionError, 'Alignment must result in less of an error.')
    end

    % Put the results together again before asigning to the surrogate
    optVect = reconstructWithFixedParameters(optVectReduced, fullProblem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display(['=== TODO_DWW: Ending parameter extraction ==='])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Si = reshapeParameters(optVect, Si, optsParE);

% Additive zero order OSM 
d = zeros(Nm,1);
if getd, d = Rfi{end} - evalSurr(xi{end},Si); end
Si.d = d;

% Additive first order OSM - ToDo


function NSMUnknowns = getNSMUnknowns()    
% NSMUnknowns - The number of SM request types represent the number of unknowns that need to be solved. 
    NSMUnknowns = 0;
    % --- getA ---
    if getA == 1
        % Diagonal
        NSMUnknowns = NSMUnknowns + 1;
    elseif getA == 2
        % Single factor
        NSMUnknowns = NSMUnknowns + 1*Nm;
    end

    % --- getB ---
    if getB == 1
        % Full
        NSMUnknowns = NSMUnknowns + Nn*Nn;
    elseif getB == 2
        % Diagonal
        NSMUnknowns = NSMUnknowns + Nn*1;
    else
        % Custom diagonal
        NSMUnknowns = NSMUnknowns + sum(getB);
    end

    % --- getc --- Additive SM
    if getc == 1
        NSMUnknowns = NSMUnknowns + Nn*1;
    end

    % --- getG --- Multiplicative ISM
    if getG == 1
        % Full
        NSMUnknowns = NSMUnknowns + Nq*Nn;
    else
        % Custom diagonal
        NSMUnknowns = NSMUnknowns + (sum(getB)*Nq);
    end

    % --- getxp --- Explicit
    if getxp == 1
        NSMUnknowns = NSMUnknowns + Nq;
    end

    % --- getd ---
    if getd == 1
        NSMUnknowns = NSMUnknowns + Nm*1;
    end

    % --- getE ---
    % Do nothing yet

    % --- getF ---
    if getF == 1
        NSMUnknowns = NSMUnknowns + 1;
    end
end % getNSMUnknowns funcation


end % buildSurr function 


% ======================================
% ========= begin subfunctions =========
% ======================================


% function [reducedProblem] = removeFixedParameters(fullProblem)
% % The full optimisation problem is reduced to only include parameters that can be changed by the optimiser.
% % The starting point for the optimisation run (x0), the left hand side inequality matrix (Aineq) as well as the
% % right hand side equality constraint (bineq) are all reduced. Lower and upper bounds are also reduced.
% % This reduces the length of optimisation runs and makes it easier to debug what is going on. The reduced
% % problem is returned in a similar form to what is passed in. Once the problem has gone through the optimiser
% % and the optimised parameters vector have been returned, the reconstructWithFixedParameters function will 
% % return the vector with the fixed parameters that were previously removed.

% % Aguments:
% %   fullProblem: The orginal problem with both modifiable and fixed parameters
% %       x0:     Initial parameter vector
% %       Aineq:  Linear inequality constraints matrix (LHS)
% %       bineq:  Linear inequality constraints vector (RHS)
% %       lb:     Lower bound vector for the optimisation parameter
% %       ub:     Upper bound vector for the optimisation parameter

% % Returns:
% %   reducedProblem: Problem for the optimiser that only has variabled that can be modified
% %       x0:     Initial parameter vector
% %       Aineq:  Linear inequality constraints matrix (LHS)
% %       bineq:  Linear inequality constraints vector (RHS)
% %       lb:     Lower bound vector for the optimisation parameter
% %       ub:     Upper bound vector for the optimisation parameter

% parameterCount = size(fullProblem.lb,1);
% assert(size(fullProblem.ub,1) == parameterCount, 'The number of parameters must all match.')
% assert(size(fullProblem.Aineq,2) == parameterCount, 'The number of parameters must all match.')
% assert(size(fullProblem.x0,1) == parameterCount, 'The number of parameters must all match.')

% reducedProblem = fullProblem;
% reducedProblem.bineq = fullProblem.bineq;

% nn = 1;
% while nn <= size(reducedProblem.lb,1)
%     if ( reducedProblem.lb(nn) == reducedProblem.ub(nn) )
%         reducedProblem.lb(nn) = [];
%         reducedProblem.ub(nn) = [];
%         % Remove any influence of the removed parameters have on the bineq (RHS vector). 
%         % The optimisation parameter affects the entire columb in the Aineq (LHS matrix).
%         reducedProblem.bineq = reducedProblem.bineq - reducedProblem.Aineq(:,nn).*reducedProblem.x0(nn);
%         reducedProblem.Aineq(:,nn) = [];
%         reducedProblem.x0(nn) = [];
%     else
%         nn = nn +1;
%     end
% end

% end % removeFixedParameters function

% ======================================

% function [optVect] = reconstructWithFixedParameters(reducedOptVect, fullProblem)
% % Builds up the optimisation vector with the fixed parameters that we removed using removeFixedParameters.
% % Fixed parameters are those whos lower and upper bounds are equal. The optimiser cannot use these values 
% % and it clutters the useful data when debugging. The original full problem is brought in for the upper and
% % lower bounds as well as the iniial values that need to be put back. 

% % Aguments:
% %   fullProblem: The orginal problem with both modifiable and fixed parameters
% %       x0:     Initial parameter vector
% %       lb:     Lower bound vector for the optimisation parameter
% %       ub:     Upper bound vector for the optimisation parameter
% %   reducedOptVect: A reduced optimisation parameter vector that only contains changeable parameters.

% % Returns:
% %   optVect: The full optimisation parameter vector with any fixed parameters filled in again.  

% parameterCount = size(fullProblem.lb,1);
% assert(size(fullProblem.ub,1) == parameterCount, 'The number of parameters must all match.')
% assert(size(fullProblem.x0,1) == parameterCount, 'The number of parameters must all match.')

% % Using zeros instead of the fullProblem.x0 for ease of debugging
% optVect = zeros(size(fullProblem.x0));

% countReduced = 1;
% for nn = 1:parameterCount
%     if ( fullProblem.lb(nn) == fullProblem.ub(nn) )
%         optVect(nn) = fullProblem.x0(nn);
%     else
%         optVect(nn) = reducedOptVect(countReduced);
%         countReduced = countReduced +1;
%     end
% end

% end % includeFixedParameters function

% ======================================

function [S] = reshapeParameters(optVect, S, optsParE)
% Reshaped the optimisation vector back into the parameters.
% Arguments:
%   optVect: optimisation vector containing all the different parameters lumped together.
%   S: The existing surrogate model
%   optsParE: options for the parameters
%       Nn;Nq;lenVect;firstPos;lastPos;errNorm;errW;
% Returns: A surrogate with its individual parameters

% Unpack the input structure
Nn = optsParE.Nn;
Nq = optsParE.Nq;
firstPos = optsParE.firstPos;
lastPos = optsParE.lastPos;
lenA = optsParE.lenVect(1);
lenB = optsParE.lenVect(2);
lenc = optsParE.lenVect(3);
lenG = optsParE.lenVect(4);
lenxp = optsParE.lenVect(5);
lenF = optsParE.lenVect(6);

% Extract individual parameters
A = diag(optVect(firstPos(1):lastPos(1)));
B = reshape(optVect(firstPos(2):lastPos(2)),sqrt(lenB),sqrt(lenB));   % Always square
c = reshape(optVect(firstPos(3):lastPos(3)),lenc,1);
G = reshape(optVect(firstPos(4):lastPos(4)),min(lenG,Nq),Nn); % Must be empty matrix if lenG == 0
xp = reshape(optVect(firstPos(5):lastPos(5)),lenxp,1);
F = reshape(optVect(firstPos(6):lastPos(6)),lenF,1); % Must be empty matrix if lenG == 0

% Update the surrogate model structure
S.A = A;
S.B = B;
S.c = c;
S.G = G;
S.xp = xp;
S.F = F;

end % reshapeParameters function

% ======================================

function e = erri(reducedOptVect,xi,Rfi,S,wk,vk,opts, fullProblem, plotFlag, plotOpts)
% Error function for optimization
% Aguments: 
%   see buildSurr for an explanation of the other arguments.
%   reducedOptVect: A reduced optimisation parameter vector that only contains changeable parameters.
%   fullProblem: The orginal problem with both modifiable and fixed parameters
%       x0:     Initial parameter vector
%       lb:     Lower bound vector for the optimisation parameter
%       ub:     Upper bound vector for the optimisation parameter
% Returns:
%   e:  The total sum of normalised error between the complex fine model and the new surrogate evaluation.
%       If multiple fine model evaluations are being used to build up a better surrogate then the error is 
%       calculated for each fine model step for each output parameter. The surrogate is evaluated for each 
%       output parameter at each of the fine model steps.

optVectn = reconstructWithFixedParameters(reducedOptVect, fullProblem);
optVect = denormaliseProblem(optVectn, fullProblem, opts)
keyboard
S = reshapeParameters(optVect, S, opts);

% Calculate the error function value
Nc = length(wk);
ec = zeros(1,Nc);
if length(opts.errW) == 1
    errW = Rfi{1}./Rfi{1};
    errW(isnan(errW)) = 1;  % In case of 0 error...
else
    errW = opts.errW;
end
diffR = {};
errorValue = {};
Rs = {};
for cc = 1:Nc
    ev = 0;
    Rs{cc} = evalSurr(xi{cc},S);
    [Nm,Np] = size(Rs{cc});
    % Errors for each output parameter (e.g. s-parameters) aggregated
    for pp = 1:Np
        diffR{cc}{pp} = errW.*(Rfi{cc}(:,pp) - Rs{cc}(:,pp));
        errorValue{cc}{pp} = norm(diffR{cc}{pp},opts.errNorm);
        ev = ev + errorValue{cc}{pp};
    end
    ec(cc) = wk(cc).*ev;
end
e = sum(ec)./Nc;

if ( plotFlag == 1 )
    plotErri(Nc, Rfi, Rs, diffR, errorValue, opts.errNorm, ec, e, plotOpts);
end

end % erri function 

% ======================================

function plotErri(Nc, Rfi, Rs, diffR, errorValue, errorNorm, ec, e, plotOpts)
% Plots the error between the fine models and the response of the surrogate.

% Parameters:
%   plotOpts:
%       plotTitle:  The graphs title.
%       yLims:  Limits for the y axis. This is used to keep the graph scaling the same as the 
%               initial alignment plot for easier comparison.
%               Format: ylim([0.05,1.5])
plotTitle = '';
if isfield(plotOpts,'plotTitle'), plotTitle = plotOpts.plotTitle; end

fig = figure();
% The number of parameters per model should not change
[Nm,Np] = size(Rs{1});
for cc = 1:Nc
    for pp = 1:Np
        diffAll = Rfi{cc}(:,pp) - Rs{cc}(:,pp);

        subplot(Nc,Np*2, (Np*(cc-1)*2) + pp*2-1), grid on, hold on
        plot(real(Rfi{cc}),'k', 'LineWidth',2)
        plot(real(Rs{cc}),'--r', 'LineWidth',2)
        plot(real(diffR{cc}{pp}),'og', 'LineWidth',2)
        plot(abs(real(diffAll)),'+m', 'LineWidth',2)
        title({[plotTitle], ...
            [' - Real:'], ...
            ['Fine model ', num2str(cc), 'of', num2str(Nc), ', Output param: ', num2str(pp), 'of', num2str(Np)], ...
            ['Outpus parameter error = ', num2str(errorValue{cc}{pp})], ...
            ['Combined normalised error = ', num2str(ec(cc))], ...
            [' using norm ', num2str(errorNorm), ', final error:', num2str(e)]})
        % legend('real(Rfi)','real(Rs)','real(diffR)','(abs(real(diffAll))')
        ylabel('Value')
        xlabel('Point')

        subplot(Nc,Np*2, (Np*(cc-1)*2) + pp*2), grid on, hold on
        plot(imag(Rfi{cc}),'k', 'LineWidth',2)
        plot(imag(Rs{cc}),'--r', 'LineWidth',2)
        plot(imag(diffR{cc}{pp}),'og', 'LineWidth',2)
        plot(abs(imag(diffAll)),'+m', 'LineWidth',2)
        title({[plotTitle], .... 
            [' - Imag:'], ...
            ['Fine model ', num2str(cc), 'of', num2str(Nc), ', Output param: ', num2str(pp), 'of', num2str(Np)], ...
            ['Outpus parameter error = ', num2str(errorValue{cc}{pp})], ...
            ['Combined normalised error = ', num2str(ec(cc))], ...
            [' using norm ', num2str(errorNorm), ', final error:', num2str(e)]})
        % legend('imag(Rfi)','imag(Rs)','imag(diffR)','(abs(imag(diffAll))')
        ylabel('Value')
        xlabel('Point')
    end
end
end % plotErri function

% ======================================

function e = erriF(Fvect,Rfi,Rc,f,opts)
% Special case error function where only F is optimized and the coarse
% model is not re-evaluated - interpolation/extrapolation is used on the
% provided coarse model response...
% This will only work for a single response at this stage

fs = Fvect(1).*f + Fvect(2);
RsComp = interp1(f,Rc,fs,'linear');
RsComp = reshape(RsComp,length(RsComp),1);
% Get rid of NaNs from shift
cleanPos = find(~isnan(RsComp));
RsCompClean = RsComp(cleanPos);
fClean = f(cleanPos);
fClean = reshape(fClean,length(fClean),1);
% Include the edge points for the interpolation
if max(fClean) < max(f) 
    RsCompClean = [RsCompClean;Rfi{1}(end)];
    fClean = [fClean;f(end)];
end
if min(fClean) > min(f)
    RsCompClean = [Rfi{1}(1);RsCompClean];
    fClean = [f(1);fClean];
end
Rs = interp1(fClean,RsCompClean,f,'spline');

diffR = Rfi{1} - reshape(Rs,length(f),1);
e = norm(diffR,opts.errNorm);

% if isequal(opts.errNorm,'L1')
%     e = sum(abs(diffR));  % Error vector [Nm,1]
% elseif isequal(opts.errNorm,'L2')
%     e = sum(abs(diffR).^2);  % Error vector [Nm,1]
% elseif isequal(opts.errNorm,'Linf')
%     e = max(abs(diffR));  % Error vector [Nm,1]
% else
%     error(['Unknown norm: ' opts.errNorm,'.  Should be L1, L2 or Linf']);
% end

end % erriF function
