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
%   Amin - minimum value for A to constrain the search space (default 0)
%   Amax - maximum value for A to constrain the search space (default inf)
%          The above values can be scalar of vectors of length Nm
%   Fmin - minimum values for F to constrain the search space [2,1] (default [0;-min(S.f)])
%   Fmax - maximum values for F to constrain the search space [2,1] (default [inf;inf])
%          Once again, equal upper and lower bounds can be used above to
%          remove a certain element from the search space
%   optsFminS - options for the fminsearch local optimizer in the PE
%   globOpt - flag to include global search (PBILreal) in the PE (Not currently working)
%   M_PBIL - number of bits per parameter in the PBIL global search (8)
%   optsPBIL - options for the PBIL global optimizer in the PE
%   errNorm - Type of error norm to use for parameter extraction (1,2,inf)
%   errW - Vector of weights (typically binary but can be any real number),
%          of length Nm, to calculate the extraction error.  Can be used to 
%          mask out regions in the response domain of less importance. 
%          Default ones. 
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

% ToDo: Implement E (first order OSM)
% ToDo: Jacobian fitting in error functions (vk)
% ToDo: Re-introduce globOpt option
% ToDo: Include more optimizer options - start with fmincon

% Preassign some variables
[lenA,lenB,lenc,lenG,lenxp,lenF] = deal(0);
[Nc,Nn,Nq,Nm] = deal(0);
[A_init,Av_init,B_init,Bv_init,c_init,G_init,Gv_init,xp_init,F_init,typeA,typeB,typec,typeG,typexp,typeF] = deal([]);
[Bmin,Bmax,cmin,cmax,Gmin,Gmax,xpmin,xpmax,fmin] = deal([]);

% Sort out formats
if iscell(xi) && iscell(Rfi) && (length(xi) == length(Rfi))   % Basic error checking - if any issue here only first point will be used
    Nc = length(xi);  % Number of input point cells
else     % Force only first point to be used, and make cell arrays
    if ~iscell(xi), xi = mat2cell(xi); end 
    if ~iscell(Rfi), Rfi = mat2cell(Rfi); end 
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
Amin = eps;
Amax = inf;
if isfield(S,'f')
    Fmin = [0,min(S.f)]';
    Fmax = [inf,inf]';
end

% Get user SM type request
if isfield(opts,'getA'), getA = opts.getA; end
if isfield(opts,'getB'), getB = opts.getB; end
if isfield(opts,'getc'), getc = opts.getc; end
if isfield(opts,'getG'), getG = opts.getG; end
if isfield(opts,'getxp'), getxp = opts.getxp; end
if isfield(opts,'getd'), getd = opts.getd; end
if isfield(opts,'getE'), getE = opts.getE; end
if isfield(opts,'getF'), getF = opts.getF; end

NSMUnknowns = getNSMUnknowns();

% Default optimization parameters
optsFminS = optimset('display','none');
globOpt = 0;
M_PBIL = 8;
optsPBIL = [];
errNorm = 2;
errW = 1;
if isfield(opts,'wk') 
    wk = opts.wk;
    if length(wk) == 1
        if (wk == 0) 
            wk = ones(1,Nc);
            % If the error function for the model fitting becomes overdetermined 
            % and can give a skewed model.
            if (Nc > NSMUnknowns)
                wk(1:end-NSMUnknowns) = 0;
            end
        else
            wk = wk.^[1:Nc];
        end
    end
else
    wk = ones(1,Nc);
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
if isfield(opts,'xpmin'), xpmin = opts.xpmin; end
if isfield(opts,'xpmax'), xpmax = opts.xpmax; end
if isfield(opts,'Amin'), Amin = opts.Amin; end
if isfield(opts,'Amax'), Amax = opts.Amax; end
if isfield(opts,'Fmin'), Fmin = opts.Fmin; end
if isfield(opts,'Fmax'), Fmax = opts.Fmax; end


% Get user optimization parameters
if isfield(opts,'optsFminS'), optsFminS = opts.optsFminS; end
% if isfield(opts,'globOpt'), globOpt = opts.globOpt; end
if isfield(opts,'globOpt'), warning('Global optimization not currently implimented - just using local constrained search'); end
if isfield(opts,'M_PBIL'), M_PBIL = opts.M_PBIL; end
if isfield(opts,'optsPBIL'), optsPBIL = opts.optsPBIL; end
if isfield(opts,'errNorm'), errNorm = opts.errNorm; end
if isfield(opts,'errW'), errW = opts.errW; end

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
lenB = Nn*Nn;
B_init = eye(Nn);
[Bv_init,Bmin,Bmax] = deal(reshape(B_init,lenB,1));
if any(getB)
    typeB = 'B';
    if isfield(S,'B')
        B_init = S.B;
    end
    % Sort out the box limits according to what is requested
    if numel(getB) == 1 && getB == 1    % Full matrix to be optimized - unconstrained in the box
        BminM = -inf.*ones(Nn,Nn);  
        BmaxM = inf.*ones(Nn,Nn);
    elseif numel(getB) == 1 && getB == 2    % Only the diagonal entries optimized - unconstrained in the box
        % Keep the rest of the entries the same as B - this allows (rarely
        % used) the off-diagonal entries to be the user supplied values
        [BminM,BmaxM] = deal(B_init);
        BminM = BminM + diag(-inf.*ones(Nn,1));
        BmaxM = BmaxM + diag(inf.*ones(Nn,1));
    elseif numel(getB) == Nn   % Only certain of the diagonal entries optimized - unconstrained in the box
        % Keep the rest of the entries the same as B - this allows (rarely
        % used) the off-diagonal entries to be the user supplied values
        [BminM,BmaxM] = deal(B_init);
        [minDiag,maxDiag] = deal(zeros(Nn,1));
        minDiag(getB == 1) = -inf;
        maxDiag(getB == 1) = inf;
        BminM = BminM + diag(minDiag);
        BmaxM = BmaxM + diag(maxDiag);
    else
        error(['Unknown getB flag: ', num2str(getB),', should be 0, 1, 2 or a bollena vector of length Nn = ', num2str(Nn)]);
    end
    Bv_init = reshape(B_init,lenB,1);
    Bmin = reshape(BminM,lenB,1);
    Bmax = reshape(BmaxM,lenB,1);
end
lenc = Nn; 
[c_init,cmin,cmax] = deal(zeros(lenc,1));
if getc
    typec = 'c';
    if isfield(S,'c')
        c_init(:,1) = S.c;
    end
    cmin = -inf.*ones(Nn,1);   % Unconstrained in the box
    cmax = inf.*ones(Nn,1);
end
% Have to provide an xp input for implicit space mapping...
if isfield(S,'xp')
    lenxp = Nq;
    [xp_init(:,1),pmin(:,1),pmax(:,1)] = deal(S.xp);    % Note name pmin/max here - xpmin reserved for linear constraints
    lenG = Nq*Nn;
    G_init = zeros(Nq,Nn); 
    [Gv_init,Gmin,Gmax] = deal(reshape(G_init,lenG,1));
    if getxp
        typexp = 'xp';
        xp_init(:,1) = S.xp;
        pmin = -inf.*ones(Nq,1);   % Unconstrained in box
        pmax = inf.*ones(Nq,1);
    end
    if any(getG)
        typeG = 'G';
        if isfield(S,'G')
            G_init = S.G;
        end
        % Sort out the box limits according to what is requested
        if numel(getG) == 1 && getG == 1    % Full matrix to be optimized - unconstrained in the box
            GminM = -inf.*ones(Nq,Nn);
            GmaxM = inf.*ones(Nq,Nn);
        elseif numel(getG) == Nn  % Only certain input parameter dependencies are kept for all implicit parameters - unconstrained in the box
            % Keep the rest of the entries the same as G
            [GminM,GmaxM] = deal(G_init);
            % Unconstrained columns where we want to allow a search
            GminM(:,getG == 1) = -inf;
            GmaxM(:,getG == 1) = inf;
        elseif all(size(getG) == [Nq,Nn])
            % Keep the rest of the entries the same as G
            [GminM,GmaxM] = deal(G_init);
            % Unconstrained entries where we want to allow a search
            GminM(getG == 1) = -inf;
            GmaxM(getG == 1) = inf;
        else
            error(['Unknown getG flag: ', num2str(getG),', should be 0, 1 or a bollean vector of length Nn = ', num2str(Nn),', or a boolean matrix of size [Nq,Nn] = [',num2str(Nq),',',num2str(Nn),']']);
        end
        Gv_init = reshape(G_init,lenG,1);
        Gmin = reshape(GminM,lenG,1);
        Gmax = reshape(GmaxM,lenG,1);
    end
else
    if getxp
        warning('Cannot perform getxp with no S.xp provided - request ignored');
    elseif any(getG)
        warning('Cannot perform getG with no S.xp provided - request ignored');
    end
    [xp_init,pmin,pmax,G_init,Gv_init,Gmin,Gmax] = deal([]);
end
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
    % Limits already handled earlier while setting up defaults and reading
    % user supplied data
    fmin = eps; % Linear constraint to force positive surrogate frequency
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
    
    if 0 && globOpt    % Only use local optimizer for FSM
        F_init = PBILreal(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),Fmin,Fmax,M_PBIL,optsPBIL);
    end
%     Fvect = fminsearchcon(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),F_init,[0, -inf],[],[-min(S.f),0;0,-1],[0;0],[],optsFminS);     % Positive multiplier, and minimum frequency
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
    Ncon = 2*Nn + 2*Nq + getF;
    LHS_mat = zeros(Ncon,length(initVect));
    
    lhsA_mat = zeros(size(A_init));        
    lhsB_mat = zeros(size(B_init));    % Basic matrix shape to use for distribution of the input vector in the LHS matrix for B
    lhsc_mat = zeros(size(c_init));    % Basic matrix shape to use for distribution of the input vector in the LHS matrix for B
    lhsG_mat = zeros(size(G_init));    % Basic matrix shape to use for distribution of the input vector in the LHS matrix for G
    lhsxp_mat = zeros(size(xp_init));
    lhsF_mat = zeros(size(F_init));
    
    RHS_vect = [-ximin;ximax;-xpmin;xpmax;-fmin];
    
    % First populate the input space limits in LHS_mat
    for xx = 1:Nn
        xA_vect = diag(lhsA_mat);       % Always zeros - no influence on input/implicit space or frequency bounds
        
        xB_mat = lhsB_mat;
        if lenB > 0;
            xB_mat(xx,:) = xi{1}';
        end
        xB_vect = reshape(xB_mat,1,lenB);
        
        xc_mat = lhsc_mat;
        if lenc > 0
            xc_mat(xx,:) = 1;
        end
        xc_vect = reshape(xc_mat,1,lenc);
        
        xG_vect = reshape(lhsG_mat,1,lenG);
        xxp_vect = reshape(lhsxp_mat,1,lenxp);
        xF_vect = reshape(lhsF_mat,1,lenF);
        
        LBrow = [xA_vect,-xB_vect,-xc_vect,xG_vect,xxp_vect,xF_vect];
        UBrow = [xA_vect,xB_vect,xc_vect,xG_vect,xxp_vect,xF_vect];
        
        LHS_mat(xx,:) = LBrow; % Lower bound row
        LHS_mat(Nn+xx,:) = UBrow; % Upper bound row
    end
    % And now the implicit parameters part
    for ii = 1:Nq
        rLB = 2*Nn + ii;    % Lower bound row
        rUB = rLB + Nq;     % Upper bound row
        
        % Always zeros
        xA_vect = diag(lhsA_mat);      
        xB_vect = reshape(lhsB_mat,1,lenB);
        xc_vect = reshape(lhsc_mat,1,lenc);
        xF_vect = reshape(lhsF_mat,1,lenF);
        
        xG_mat = lhsG_mat;
        if lenG > 0
            xG_mat(ii,:) = xi{1}';
        end
        xG_vect = reshape(xG_mat,1,lenG);
        
        xxp_mat = lhsxp_mat;
        if lenxp > 0
            xxp_mat(ii,:) = 1;
        end
        xxp_vect = reshape(xxp_mat,1,lenxp);
                
        LBrow = [xA_vect,xB_vect,xc_vect,-xG_vect,-xxp_vect,xF_vect];
        UBrow = [xA_vect,xB_vect,xc_vect,xG_vect,xxp_vect,xF_vect];
        
        LHS_mat(rLB,:) = LBrow; % Lower bound row
        LHS_mat(rUB,:) = UBrow; % Upper bound row
    end
    % And frequency shifts
    if getF
        LHS_mat(2*(Nn+Nq)+1,end-1:end) = [-min(S.f),-1];
    end
    
    if 0 && globOpt     % Not implimented yet
        % Start with global search to get initial value
        [initVect] = PBILreal(@(optVect) erri(optVect,xi,Rfi,S,wk,vk,optsParE),minVect,maxVect,M_PBIL,optsPBIL);
    end
%     optVect = fminsearch(@(optVect) erri(optVect,xi,Rfi,S,wk,vk,optsParE),initVect,optsFminS);
    optVect = fminsearchcon(@(optVect) erri(optVect,xi,Rfi,S,wk,vk,optsParE),initVect,minVect,maxVect,LHS_mat,RHS_vect,[],optsFminS);
end
% Rebuild the individual parameters from the vector
A = diag(optVect(firstPos(1):lastPos(1)));
B = reshape(optVect(firstPos(2):lastPos(2)),sqrt(lenB),sqrt(lenB));   % Always square
c = reshape(optVect(firstPos(3):lastPos(3)),lenc,1);
G = reshape(optVect(firstPos(4):lastPos(4)),min(lenG,Nq),Nn); % Must be empty matrix if lenG == 0
xp = reshape(optVect(firstPos(5):lastPos(5)),lenxp,1);
F = reshape(optVect(firstPos(6):lastPos(6)),lenF,1); % Must be empty matrix if lenG == 0

Si.A = A; 
Si.B = B; 
Si.c = c; 
Si.G = G; 
Si.xp = xp; 
Si.F = F; 

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
    if getd == 1
        NSMUnknowns = NSMUnknowns + 1;
    end
end % end getNSMUnknowns

end % end buildSurr


% Error function for optimization
function e = erri(optVect,xi,Rfi,S,wk,vk,opts)

% Unpack the input structure
Nn = opts.Nn;
Nq = opts.Nq;
firstPos = opts.firstPos;
lastPos = opts.lastPos;
lenA = opts.lenVect(1);
lenB = opts.lenVect(2);
lenc = opts.lenVect(3);
lenG = opts.lenVect(4);
lenxp = opts.lenVect(5);
lenF = opts.lenVect(6);

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

% Calculate the error function value
Nc = length(wk);
ec = zeros(1,Nc);
if length(opts.errW) == 1
    errW = Rfi{1}./Rfi{1};
    errW(isnan(errW)) = 1;  % In case of 0 error...
else
    errW = opts.errW;
end
% Count over points
for cc = 1:Nc
    Rs = evalSurr(xi{cc},S);
    [Nm,Np] = size(Rs);
    ev = 0;
    % Errors for each of the output parameters aggregated
    for pp = 1:Np
        diffR = errW.*(Rfi{cc}(:,pp) - Rs(:,pp));
        ev = ev + norm(diffR,opts.errNorm);
    end
    ec(cc) = wk(cc).*ev;
end
e = sum(ec)./Nc;
end


% Special case error function where only F is optimized and the coarse
% model is not re-evaluated - interpolation/extrapolation is used on the
% provided coarse model response...
% This will only work for a single response at this stage
function e = erriF(Fvect,Rfi,Rc,f,opts)
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

end

