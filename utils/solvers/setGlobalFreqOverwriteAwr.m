function setGlobalFreqOverwriteAwr(overwiteFlag)
% This is a workaround for a bug in AWR. The frequency items cannot be queried and set 
% correctly. If that functionality becomes available then this global hack can be removed.
% Clearing and rewriting all the frequencies in the AWR model is incredibly slow and this
% should only be used if it is really needed with frequency space-mapping.
% https://uk.mathworks.com/help/matlab/ref/global.html

% Argument: a flag to set the global overwrite flag for AWR.

global globalFreqOverwriteAwr
globalFreqOverwriteAwr = overwiteFlag;

end % --- setGlobalFreqOverwriteAwr