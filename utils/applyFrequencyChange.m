function [Rs] = applyFrequencyChange(origninalFreq, Fvect, Rc)
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

fs = Fvect(1).*origninalFreq + Fvect(2);

RsComp = interp1(origninalFreq, Rc, fs, 'spline');
RsComp = reshape(RsComp, length(RsComp),1);
% Get rid of NaNs from shift
cleanPos = find(~isnan(RsComp));
RsCompClean = RsComp(cleanPos);
fClean = origninalFreq(cleanPos);
fClean = reshape(fClean, length(fClean), 1);
% Include the edge points for the interpolation
if max(fClean) < max(origninalFreq)
    RsCompClean = [RsCompClean; Rc(end)];
    fClean = [fClean; origninalFreq(end)];
end
if min(fClean) > min(origninalFreq)
    RsCompClean = [Rc(1); RsCompClean];
    fClean = [origninalFreq(1); fClean];
end
Rs = interp1(fClean, RsCompClean, origninalFreq, 'spline');

end % applyFrequencyChange function