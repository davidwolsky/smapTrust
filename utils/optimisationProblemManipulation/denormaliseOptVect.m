function x0 = denormaliseOptVect(optVectn, originalProb, optsParE)

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
x0_An  = optVectn(firstPos(1):lastPos(1));
x0_Bn  = optVectn(firstPos(2):lastPos(2));
x0_cn  = optVectn(firstPos(3):lastPos(3));
x0_Gn  = optVectn(firstPos(4):lastPos(4));
x0_xpn = optVectn(firstPos(5):lastPos(5));
x0_Fn  = optVectn(firstPos(6):lastPos(6));

% --- Extract bounds ---

xmin = reshape(-originalProb.bineq(1   : Nn),   Nn,1);
xmax = reshape(originalProb.bineq(Nn+1 : 2*Nn), Nn,1);

xpmin = reshape(-originalProb.bineq(2*Nn+1   : 2*Nn+Nq),   Nq,1);
xpmax = reshape(originalProb.bineq(2*Nn+Nq+1 : 2*Nn+2*Nq), Nq,1);

fmin = -originalProb.bineq(end-1);
fmax = originalProb.bineq(end);

deltax = xmax - xmin;
deltaxp = xpmax - xpmin;

deltaf = fmax - fmin;

% --- Denormalise x0 --- 
x0_A = x0_An;

x0_B = x0_Bn;

x0_c = x0_cn.*deltax + xmin;

GnFact = 10;
x0_G = x0_Gn/GnFact;

x0_xp = x0_xpn.*deltaxp + xpmin;

F1 = x0_Fn(1);
F2 = x0_Fn(2).*deltaf + fmin;
x0_F = [F1; F2];

x0 = [x0_A(:); x0_B(:); x0_c(:); x0_G(:); x0_xp(:); x0_F(:)];

end % denormaliseOptVect function