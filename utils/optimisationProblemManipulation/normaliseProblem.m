function nProb = normaliseProblem(prob, optsParE)

% Normalises before parameter extraction takes place. SM parameters that are linked to model variabels
% are normalised but multiplicative parameters are left as is or scaled slightly. 
% Arguments:
%   optsParE:   Options for parameter extraction.
%       - Nn:   Number of design variabled.
%       - Nq:   Number of implicit variables.          
%       - lenVect:  Vector containing position/lengths of each SM parameter.
%       - firstPos: First position of each SM parameter set. 
%       - lastPos:  Last position of each SM parameter set. 
%       - errNorm:  Error norm in use.
%       - errW: Vector of weights (typically binary but can be any real number),
%               of length Nm, to calculate the extraction error.  Can be used to 
%               mask out regions in the response domain of less importance. 
%               Default ones.
% Returns:
%   nProb:  Normalised problem: x0, ub, lb, Aineq and bineq.

parameterCount = size(prob.lb,1);
assert(size(prob.ub,1) == parameterCount, 'The number of parameters must all match.')
assert(size(prob.x0,1) == parameterCount, 'The number of parameters must all match.')
nProb = prob;

Nn = optsParE.Nn;
Nq = optsParE.Nq;
firstPos = optsParE.firstPos;
lastPos = optsParE.lastPos;
% lenA = optsParE.lenVect(1);
% lenB = optsParE.lenVect(2);
% lenc = optsParE.lenVect(3);
% lenG = optsParE.lenVect(4);
% lenxp = optsParE.lenVect(5);
% lenF = optsParE.lenVect(6);

% Extract individual parameters and bounds
x0_A  = prob.x0(firstPos(1):lastPos(1));
x0_B  = prob.x0(firstPos(2):lastPos(2));
x0_c  = prob.x0(firstPos(3):lastPos(3));
x0_G  = prob.x0(firstPos(4):lastPos(4));
x0_xp = prob.x0(firstPos(5):lastPos(5));
x0_F  = prob.x0(firstPos(6):lastPos(6));

x0_Amin = prob.lb(firstPos(1):lastPos(1));
x0_Bmin = prob.lb(firstPos(2):lastPos(2));
x0_cmin = prob.lb(firstPos(3):lastPos(3));
x0_Gmin = prob.lb(firstPos(4):lastPos(4));
x0_pmin = prob.lb(firstPos(5):lastPos(5));
x0_Fmin = prob.lb(firstPos(6):lastPos(6));

x0_Amax = prob.ub(firstPos(1):lastPos(1));
x0_Bmax = prob.ub(firstPos(2):lastPos(2));
x0_cmax = prob.ub(firstPos(3):lastPos(3));
x0_Gmax = prob.ub(firstPos(4):lastPos(4));
x0_pmax = prob.ub(firstPos(5):lastPos(5));
x0_Fmax = prob.ub(firstPos(6):lastPos(6));

% --- Extract bounds ---

% RHS_vect = [-ximin; ximax; -xpmin; xpmax; -fmin; fmax];
xmin = reshape(-prob.bineq(1   : Nn),   Nn,1);
xmax = reshape(prob.bineq(Nn+1 : 2*Nn), Nn,1);

xpmin = reshape(-prob.bineq(2*Nn+1   : 2*Nn+Nq),   Nq,1);
xpmax = reshape(prob.bineq(2*Nn+Nq+1 : 2*Nn+2*Nq), Nq,1);

fmin = -prob.bineq(2*Nn+2*Nq+1);
fmax = prob.bineq(2*Nn+2*Nq+2);

deltax = xmax - xmin;
deltaxp = xpmax - xpmin;

deltaf = fmax - fmin;

% --- Normalise x0 and bounds --- 
x0_An = x0_A;
x0_Albn = x0_Amin;
x0_Aubn = x0_Amax;

x0_Bn = x0_B;
x0_Blbn = x0_Bmin;
x0_Bubn = x0_Bmax;

x0_cn   = (x0_c - xmin)    ./ (deltax);
x0_clbn = (x0_cmin - xmin) ./ (deltax);
x0_cubn = (x0_cmax - xmin) ./ (deltax);

GnFact = 10;
x0_Gn   = GnFact*x0_G;
x0_Glbn = GnFact*x0_Gmin;
x0_Gubn = GnFact*x0_Gmax;

x0_xpn  = (x0_xp - xpmin)  ./ (deltaxp);
x0_plbn = (x0_pmin- xpmin) ./ (deltaxp);
x0_pubn = (x0_pmax- xpmin) ./ (deltaxp);

x0_F1n = x0_F(1);
x0_F1lbn = x0_Fmin(1);
x0_F1ubn = x0_Fmax(1);

x0_F2n   = (x0_F(2) - fmin)    ./ (deltaf);
x0_F2lbn = (x0_Fmin(2) - fmin) ./ (deltaf);
x0_F2ubn = (x0_Fmax(2) - fmin) ./ (deltaf);

x0_Fn = [x0_F1n; x0_F2n];
x0_Flbn = [x0_F1lbn; x0_F2lbn];
x0_Fubn = [x0_F1ubn; x0_F2ubn];

nProb.x0 = [x0_An(:); x0_Bn(:); x0_cn(:); x0_Gn(:); x0_xpn(:); x0_Fn(:)];
nProb.lb = [x0_Albn; x0_Blbn; x0_clbn; x0_Glbn; x0_plbn; x0_Flbn];
nProb.ub = [x0_Aubn; x0_Bubn; x0_cubn; x0_Gubn; x0_pubn; x0_Fubn];

% --- Normalise inequalities ---
nProb.Aineq = prob.Aineq;

deltax_mat = diag(deltax);
% Item 3 is c
nProb.Aineq(1:Nn,      firstPos(3):lastPos(3)) = -deltax_mat;
nProb.Aineq(Nn+1:2*Nn, firstPos(3):lastPos(3)) = deltax_mat;

deltaxp_mat = diag(deltaxp);
% Item 5 is xp
nProb.Aineq(2*Nn+1:   2*Nn+Nq,   firstPos(5):lastPos(5)) = -deltaxp_mat;
nProb.Aineq(2*Nn+Nq+1:2*Nn+2*Nq, firstPos(5):lastPos(5)) = deltaxp_mat;

% Last is F2
nProb.Aineq(end-1,  end) = -deltaf;
nProb.Aineq(end,    end) = deltaf;

RHSnTerm = [-xmin; xmin; -xpmin; xpmin; -fmin; fmin];
nProb.bineq = prob.bineq - RHSnTerm;

end % normaliseProblem main function