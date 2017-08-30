% TODO_DWW: Remove
function problem = denormaliseProblem(nProb, originalProb, optsParE)

parameterCount = size(nProb.lb,1);
assert(size(nProb.ub,1) == parameterCount, 'The number of parameters must all match.')
assert(size(nProb.x0,1) == parameterCount, 'The number of parameters must all match.')

problem = nProb;

firstPos = optsParE.firstPos;
lastPos = optsParE.lastPos;

keyboard
% --- A - firstPos(1) ---

% --- B - firstPos(2) ---

% --- c - firstPos(3) ---
denormaliseInBounds(firstPos(3), lastPos(3))

% --- G - firstPos(4) ---

% --- xp - firstPos(5) ---

% --- F - firstPos(6) ---


% nProblem.lb = prob.lb - prob.lb
% nProblem.ub = prob.ub ./ prob.ub
% nProblem.x0 = (prob.x0 - prob.lb)./(prob.ub - prob.lb)

% ======================================
% ========= begin subfunctions =========
% ======================================

function denormaliseInBounds(startValue, endValue)
for nn = startValue:endValue
    if ( isfinite(originalProb.ub(nn)) )
        problem.lb(nn) = originalProb.lb(nn);
        problem.ub(nn) = originalProb.ub(nn);
        problem.x0(nn) = (nProb.x0(nn))*(originalProb.ub(nn) - originalProb.lb(nn)) + originalProb.lb(nn);
    end
end
end % denormaliseInBounds function

end % denormaliseProblem function