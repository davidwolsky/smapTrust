function [optVect, fval, exitflag, output] = doOptimisation(problem)

switch problem.solver
    % ----- Local solvers ----- 
    case 'fmincon'
        [optVect, fval, exitflag, output] = fmincon(problem);
    case 'fminsearchcon'
        problem.options.OutputFcn = [];
        [optVect, fval, exitflag, output] = fminsearchcon(problem.objective, problem.x0, ...
                                                          problem.lb, problem.ub, ...
                                                          problem.Aeq, problem.beq, ...
                                                          problem.nonlcon, problem.options);
    case 'fminsearch'
        [optVect, fval, exitflag, output] = fminsearch(problem);
    % ----- Global solvers ----- 
    case 'ga'
        [optVect, fval, exitflag, output] = ga(problem);
    case 'patternsearch'
        [optVect, fval, exitflag, output] = patternsearch(problem);
    otherwise
        error(['Unknows solver type (', problem.solver,')'])
end

end % doOptimisation function





