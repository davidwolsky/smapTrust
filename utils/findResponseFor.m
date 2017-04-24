function Ri = findResponseFor(R, goalResType)

% Function to find a corresponding response for the goal result type in question (goalResType).
% R is a cell array of structures containing the response in R.r, the type R.t, and the
% (optional) domain (typically frequency) in R.f.
% R can also be a structure if only one type of response is considered.
% If no response type or domain is specified, R may also be a vector.
% A Goal can contain (typically a subset of OPTopts used in the main function):
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'S11_complex'
%               'S11_real'
%               'S11_imag'
%               'S11_dB'
%               'S11_abs'
%               'S11_angle'
%               'S11_deg' - degrees
%               'Gen' - general 

foundMathingType = false;
Nr = length(R);
tt = 1;
while tt <= Nr
    if isfield(R{tt},'t') && strncmp(R{tt}.t, goalResType, 3)
        Ri = convertResponse(R{tt}, goalResType);
        foundMathingType = true;
        break;
    end
    tt = tt + 1;
end
assert(foundMathingType, ['No matching result type was found for the specified goalResType.', goalResType, '.  R{:}.t = ', R{:}.t])

end % function findResponseFor