function nProblem = normaliseProblem(prob, optsParE)

parameterCount = size(prob.lb,1);
assert(size(prob.ub,1) == parameterCount, 'The number of parameters must all match.')
assert(size(prob.x0,1) == parameterCount, 'The number of parameters must all match.')
nProblem = prob;

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
A  = [prob.x0(firstPos(1):lastPos(1))];
B  = [prob.x0(firstPos(2):lastPos(2))];
c  = [prob.x0(firstPos(3):lastPos(3))];
G  = [prob.x0(firstPos(4):lastPos(4))];
xp = [prob.x0(firstPos(5):lastPos(5))];
F  = [prob.x0(firstPos(6):lastPos(6))];

Amin = [prob.lb(firstPos(1):lastPos(1))];
Bmin = [prob.lb(firstPos(2):lastPos(2))];
cmin = [prob.lb(firstPos(3):lastPos(3))];
Gmin = [prob.lb(firstPos(4):lastPos(4))];
pmin = [prob.lb(firstPos(5):lastPos(5))];
Fmin = [prob.lb(firstPos(6):lastPos(6))];

Amax = [prob.ub(firstPos(1):lastPos(1))];
Bmax = [prob.ub(firstPos(2):lastPos(2))];
cmax = [prob.ub(firstPos(3):lastPos(3))];
Gmax = [prob.ub(firstPos(4):lastPos(4))];
pmax = [prob.ub(firstPos(5):lastPos(5))];
Fmax = [prob.ub(firstPos(6):lastPos(6))];

% --- Extract bounds ---
xmin(:) = -prob.bineq(1 : Nn);
xmax(:) = prob.bineq(Nn+1 : 2*Nn);

xpmin(:) = -prob.bineq(2*Nn+1 : 2*Nn+Nq);
xpmax(:) = prob.bineq(2*Nn+Nq+1 : 2*Nn+2*Nq);

keyboard
% TODO_DWW: Find a better way to turn this off
if length(prob.bineq) > 2*Nn+2*Nq
    fmin = -prob.bineq(2*Nn+2*Nq+1);
    fmax = prob.bineq(2*Nn+2*Nq+2);
end

deltax = xmax - xmin;
deltaxp = xpmax - xpmin;

% TODO_DWW: Find a better way to turn this off
if length(prob.bineq) > 2*Nn+2*Nq
    deltaf = fmax - fmin;
% else 
    % deltaf = 0
end

% --- Normalise x0 and bounds --- 
An = A;
Albn = Amin;
Aubn = Amax;

Bn = B;
Blbn = Bmin;
Bubn = Bmax;

cn = (c - xmin) ./ (deltax);
clbn = (cmin- xmin) ./ (deltax);
cubn = (cmax- xmin) ./ (deltax);

GnFact = 10;
Gn = GnFact*G;
Glbn = GnFact*Gmin;
Gubn = GnFact*Gmax;

xpn = (xp - xpmin) ./ (deltaxp);
plbn = (pmin- xpmin) ./ (deltaxp);
pubn = (pmax- xpmin) ./ (deltaxp);


% TODO_DWW: Find a better way to turn this off
if length(prob.bineq) > 2*Nn+2*Nq
    F1n = F(1);
    F1lbn = Fmin(1);
    F1ubn = Fmax(1);

    F2n = (F(2) - fmin) ./ (deltaf);
    F2lbn = (Fmin(2) - fmin) ./ (deltaf);
    F2ubn = (Fmax(2) - fmin) ./ (deltaf);

    Fn = [F1n;F2n];
    Flbn = [F1lbn;F2lbn];
    Fubn = [F1ubn;F2ubn];
else
    Fn = [1;0]
    Flbn = [1;0]
    Fubn = [1;0]
end

nProb.x0 = [An(:); Bn(:); cn(:); Gn(:); xpn(:); Fn(:)];
nProb.lb = [Albn; Blbn; clbn; Glbn; plbn; Flbn];
nProb.ub = [Aubn; Bubn; cubn; Gubn; pubn; Fubn];

% --- Normalise inequalities ---
nProb.Aineq = prob.Aineq;

deltax_mat = diag(deltax);
nProb.Aineq(1:Nn,      firstPos(3):lastPos(3)) = -deltax_mat;
nProb.Aineq(Nn+1:2*Nn, firstPos(3):lastPos(3)) = deltax_mat;

deltaxp_mat = diag(deltaxp);
rxpStart = 2*Nn+1;
nProb.Aineq(rxpStart:rxpStart+Nq,        firstPos(5):lastPos(5)) = -deltaxp_mat;
nProb.Aineq(rxpStart+Nq+1:rxpStart+2*Nq, firstPos(5):lastPos(5)) = deltaxp_mat;


% TODO_DWW: Find a better way to turn this off
if length(prob.bineq) > 2*Nn+2*Nq
    nProb.Aineq(end-1,  end) = -deltaf;
    nProb.Aineq(end,    end) = deltaf;
else
    nProb.Aineq(end-1,  end) = 0;
    nProb.Aineq(end,    end) = 0;
    fmin = []
end

keyboard
RHSnTerm = [-xmin; xmin; -xpmin; xpmin; -fmin; fmin];
nProb.bineq = prob.bineq - RHSnTerm;

end % normaliseProblem main function