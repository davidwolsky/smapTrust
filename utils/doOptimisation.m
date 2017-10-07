function [optVect, fval, exitflag, output] = doOptimisation(problem)

switch problem.solver
    % ----- Local solvers ----- 
    case 'fmincon'
        [optVect, fval, exitflag, output] = fmincon(problem);
    % ----- Global solvers ----- 
    case 'ga'
        [optVect, fval, exitflag, output] = ga(problem);
    otherwise
        error(['Unknows solver type (', problem.solver,')'])
end

end % doOptimisation function