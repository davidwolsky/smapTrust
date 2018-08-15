function [R]=interpolateInf(R, f)
% TODO_DWW: Explain this...
% TODO_DWW: Maybe add a NaN one to simplify applyFreqChange function... 
% Get rid of inf from data/dB.

if any(isinf(R))
%     keyboard
    infPos = find(isinf(R));
    cleanPos = find(~isinf(R));
    RClean = R(cleanPos);
    fClean = f(cleanPos);
    infInterp = griddedInterpolant(fClean, RClean);
    R(infPos) = infInterp(f(infPos));
    assert(length(R)==length(f), 'Vectors must match size.')
end

end % interpolateInf function