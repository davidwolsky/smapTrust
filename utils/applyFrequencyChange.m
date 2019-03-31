function [Rs] = applyFrequencyChange(origininalFreq, Fvect, Rc)
% This applies a scaling and offset to the frequency and then interpiplates the response to match the 
% new frequency axis. The new responce is returned.

% Returns:
%   Rs: The new (surrogate) response after it has been adapted.

% Inputs:
%   origninalFreq:  The inital vector of frequencies.    
%   Fvect:  A vector housing the scaling and offset to be applied to the frequency. First item scaling, second offset.
%           Format -> fs = Fvect(1).*origninalFreq + Fvect(2);
%           These are typically the optimised parameters.
%   Rc: Coarse model response. This is interpolated on the new axis. Also used for filling in end values.
%       Assuming that this does not include phase (pass in absolute value of Rc)

assert(~any(isinf(Rc)), 'Need to have clean data to work with. Suggest interpolateInf method.')

fs = Fvect(1).*origininalFreq + Fvect(2);
% keyboard
% TODO_DWW: This has to be linear.
% TODO_DWW: Comment about NaN and extrapolation
RsComp = interp1(origininalFreq, Rc, fs, 'linear');
RsComp = reshape(RsComp, length(RsComp),1);

% Get rid of NaNs from shift
cleanPos = find(~isnan(RsComp));
assert(~isempty(cleanPos), 'Cannot work with only NaN values.')
RsCompClean = RsComp(cleanPos);
fClean = origininalFreq(cleanPos);
fClean = reshape(fClean, length(fClean), 1);

% Include the edge points for the interpolation
if max(fClean) < max(origininalFreq)
    RsCompClean = [RsCompClean; Rc(end)];
    fClean = [fClean; origininalFreq(end)];
end
if min(fClean) > min(origininalFreq)
    RsCompClean = [Rc(1); RsCompClean];
    fClean = [origininalFreq(1); fClean];
end

% TODO_DWW: Comment about NaN and extrapolation and shape fitting
Rs = interp1(fClean, RsCompClean, origininalFreq, 'pchip');

end % applyFrequencyChange function