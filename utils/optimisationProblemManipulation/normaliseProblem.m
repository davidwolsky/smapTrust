% TODO_DWW: Remove
function nProblem = normaliseProblem(prob, optsParE)

parameterCount = size(prob.lb,1);
assert(size(prob.ub,1) == parameterCount, 'The number of parameters must all match.')
assert(size(prob.x0,1) == parameterCount, 'The number of parameters must all match.')
nProblem = prob;

firstPos = optsParE.firstPos;
lastPos = optsParE.lastPos;

% keyboard

% --- A - firstPos(1) ---

% --- B - firstPos(2) ---

% --- c - firstPos(3) ---
normaliseInBounds(firstPos(3), lastPos(3))

% --- G - firstPos(4) ---

% --- xp - firstPos(5) ---

% --- F - firstPos(6) ---



% nProblem.lb = prob.lb - prob.lb
% nProblem.ub = prob.ub ./ prob.ub
% nProblem.x0 = (prob.x0 - prob.lb)./(prob.ub - prob.lb)

% ======================================
% ========= begin subfunctions =========
% ======================================

function normaliseInBounds(startValue, endValue)
for nn = startValue:endValue
    if ( isfinite(prob.ub(nn)) )
        nProblem.lb(nn) = prob.lb(nn) - prob.lb(nn);
        nProblem.ub(nn) = prob.ub(nn) ./ prob.ub(nn);
        nProblem.x0(nn) = (prob.x0(nn) - prob.lb(nn))./(prob.ub(nn) - prob.lb(nn));

        % TODO_DWW: 
        % Aineq
        % bineq
    end
end
end % normaliseInBounds function

end % normaliseProblem main function