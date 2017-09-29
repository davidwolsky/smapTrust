function [] = createLog(filename, space)

% Creates a log file from the variables stored in space.
%   filename:   The filename for the log file that stores the workspace.
%   space:      A grouping variable that houses all the values that should be stored in the log file.

save(filename, 'space')

end % createLog function