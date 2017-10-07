function Ri = findResponseFor(R, goalResType)

% Function to find a corresponding response for the goal result type in question (goalResType).
% R is a cell array of structures containing the response in R.r, the type R.t, and the
% (optional) domain (typically frequency) in R.f.
% R can also be a structure if only one type of response is considered.
% If no response type or domain is specified, R may also be a vector.
% A Goal can contain (typically a subset of OPTopts used in the main function):
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'Sb,a_complex'
%               'Sb,a_real'
%               'Sb,a_imag'
%               'Sb,a_dB'
%               'Sb,a_abs'
%               'Sb,a_angle'
%               'Sb,a_deg' - degrees
%               'Gen' - general 

foundMatchingType = false;
Nr = length(R);
tt = 1;
while tt <= Nr
    if isfield(R{tt},'t') && ( strcmp(R{tt}.t, goalResType) || strncmp(R{tt}.t, goalResType, find(goalResType=='_')-1) )
        Ri = convertResponse(R{tt}, goalResType);
        foundMatchingType = true;
        break;
    end
    tt = tt + 1;
end
% keyboard
% assert(foundMatchingType, ['No matching result type was found for the specified goalResType.', goalResType, '.  R{:}.t = ', R{:}.t])
assert(foundMatchingType, ['No matching result type was found for the specified goalResType: ', goalResType])

end % function findResponseFor