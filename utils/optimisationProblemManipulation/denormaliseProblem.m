% TODO_DWW: Rename denormaliseOptVect?
function x0 = denormaliseProblem(optVectn, originalProb, optsParE)

% parameterCount = size(nProb.lb,1);
% assert(size(nProb.ub,1) == parameterCount, 'The number of parameters must all match.')
% assert(size(nProb.x0,1) == parameterCount, 'The number of parameters must all match.')

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
An  = [optVectn(firstPos(1):lastPos(1))];
Bn  = [optVectn(firstPos(2):lastPos(2))];
cn  = [optVectn(firstPos(3):lastPos(3))];
Gn  = [optVectn(firstPos(4):lastPos(4))];
xpn = [optVectn(firstPos(5):lastPos(5))];
Fn  = [optVectn(firstPos(6):lastPos(6))];

% --- Extract bounds ---
xmin(:) = -originalProb.bineq(1 : Nn);
xmax(:) = originalProb.bineq(Nn+1 : 2*Nn);

xpmin(:) = -originalProb.bineq(2*Nn+1 : 2*Nn+Nq);
xpmax(:) = originalProb.bineq(2*Nn+Nq+1 : 2*Nn+2*Nq);

% TODO_DWW: Find a better way to turn this off
if length(originalProb.bineq) > 2*Nn+2*Nq
    fmin = -originalProb.bineq(end-1);
    fmax = originalProb.bineq(end);
end

deltax = xmax - xmin;
deltaxp = xpmax - xpmin;

% TODO_DWW: Find a better way to turn this off
if length(originalProb.bineq) > 2*Nn+2*Nq
    deltaf = fmax - fmin;
end

% --- Denormalise x0 --- 
A = An;

B = Bn;

c = cn.*deltax + xmin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                c is fucked up                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GnFact = 10;
G = Gn/GnFact;

xp = xpn.*deltaxp + xpmin;

keyboard
% TODO_DWW: Find a better way to turn this off
if length(originalProb.bineq) > 2*Nn+2*Nq
    F1 = Fn(1);
    F2 = Fn(2).*deltaf + fmin;
    F = [F1;F2];
else
    F = [1;0];
end

x0 = [A(:); B(:); c(:); G(:); xp(:); F(:)];

end % denormaliseProblem function