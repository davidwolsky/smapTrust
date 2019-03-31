function [reducedProblem] = removeFixedParameters(fullProblem)
% The full optimisation problem is reduced to only include parameters that can be changed by the optimiser.
% The starting point for the optimisation run (x0), the left hand side inequality matrix (Aineq) as well as the
% right hand side equality constraint (bineq) are all reduced. Lower and upper bounds are also reduced.
% This reduces the length of optimisation runs and makes it easier to debug what is going on. The reduced
% problem is returned in a similar form to what is passed in. Once the problem has gone through the optimiser
% and the optimised parameters vector have been returned, the reconstructWithFixedParameters function will 
% return the vector with the fixed parameters that were previously removed.

% Aguments:
%   fullProblem: The orginal problem with both modifiable and fixed parameters
%       x0:     Initial parameter vector
%       Aineq:  Linear inequality constraints matrix (LHS)
%       bineq:  Linear inequality constraints vector (RHS)
%       lb:     Lower bound vector for the optimisation parameter
%       ub:     Upper bound vector for the optimisation parameter

% Returns:
%   reducedProblem: Problem for the optimiser that only has variabled that can be modified
%       x0:     Initial parameter vector
%       Aineq:  Linear inequality constraints matrix (LHS)
%       bineq:  Linear inequality constraints vector (RHS)
%       lb:     Lower bound vector for the optimisation parameter
%       ub:     Upper bound vector for the optimisation parameter

parameterCount = size(fullProblem.lb,1);
assert(size(fullProblem.ub,1) == parameterCount, 'The number of parameters must all match.')
assert(size(fullProblem.Aineq,2) == parameterCount, 'The number of parameters must all match.')
assert(size(fullProblem.x0,1) == parameterCount, 'The number of parameters must all match.')

reducedProblem = fullProblem;
reducedProblem.bineq = fullProblem.bineq;

nn = 1;
while nn <= size(reducedProblem.lb, 1)
    if ( reducedProblem.lb(nn) == reducedProblem.ub(nn) )
        reducedProblem.lb(nn) = [];
        reducedProblem.ub(nn) = [];
        % Remove any influence of the removed parameters have on the bineq (RHS vector). 
        % The optimisation parameter affects the entire columb in the Aineq (LHS matrix).
        % keyboard
        reducedProblem.bineq = reducedProblem.bineq - reducedProblem.Aineq(:,nn).*reducedProblem.x0(nn);
        reducedProblem.Aineq(:,nn) = [];
        reducedProblem.x0(nn) = [];
    else
        nn = nn +1;
    end
end

end % removeFixedParameters function
