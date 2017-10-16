function [space] = getPrepopulatedSpaceFrom(filename)

% Load the 'space' varaible from a file. This is a grouping of variables of interest from a previous 
% workspace/set of runs. This is typically used for loading the Ti.xi_all, Ti.Rfi_all, Ti.costF_all 
% variables from the previous space. They are used to build a better surrogate.
%   filename: the filename used to load variables from.
% Returns 
%   space: a grouping of variables loaded from a log file.

load(filename, 'space')

end % getPrepopulatedSpaceFrom function
