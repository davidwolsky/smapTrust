function [optVect] = reconstructWithFixedParameters(reducedOptVect, fullProblem)
% Builds up the optimisation vector with the fixed parameters that we removed using removeFixedParameters.
% Fixed parameters are those whos lower and upper bounds are equal. The optimiser cannot use these values 
% and it clutters the useful data when debugging. The original full problem is brought in for the upper and
% lower bounds as well as the iniial values that need to be put back. 

% Aguments:
%   fullProblem: The orginal problem with both modifiable and fixed parameters
%       x0:     Initial parameter vector
%       lb:     Lower bound vector for the optimisation parameter
%       ub:     Upper bound vector for the optimisation parameter
%   reducedOptVect: A reduced optimisation parameter vector that only contains changeable parameters.

% Returns:
%   optVect: The full optimisation parameter vector with any fixed parameters filled in again.  

parameterCount = size(fullProblem.lb,1);
assert(size(fullProblem.ub,1) == parameterCount, 'The number of parameters must all match.')
assert(size(fullProblem.x0,1) == parameterCount, 'The number of parameters must all match.')

% Using zeros instead of the fullProblem.x0 for ease of debugging
optVect = zeros(size(fullProblem.x0));

countReduced = 1;
for nn = 1:parameterCount
    if ( fullProblem.lb(nn) == fullProblem.ub(nn) )
        optVect(nn) = fullProblem.x0(nn);
    else
        optVect(nn) = reducedOptVect(countReduced);
        countReduced = countReduced +1;
    end
end

end % includeFixedParameters function