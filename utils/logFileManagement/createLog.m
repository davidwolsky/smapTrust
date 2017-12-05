function [] = createLog(filename, space)

% Creates a log file from the variables stored in space.
%   filename:   The filename for the log file that stores the workspace. The current 
%               time and date will be appended.
%               The file is saved from where the script is called from. 
%   space:      A grouping variable that houses all the values that should be stored in the log file.

save([filename, datestr(datetime('now'), '_yyyymmdd_HHMM'), '.mat'], 'space')

end % createLog function