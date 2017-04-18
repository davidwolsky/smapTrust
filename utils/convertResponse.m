function Ri = convertResponse(R, goalResType)

% Function to convert a response (R) from one type (R.t) to a specified goalResType. 
% R is a cell array of structures containing the response in R.r, the type R.t, and the
% (optional) domain (typically frequency) in R.f.
% R can also be a structure if only one type of response is considered.
% If no response type or domain is specified, R may also be a vector.
% A Goal can contain (typically a subset of OPTopts used in the main function):
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
% TODO_DWW: Update comments everywhere to have the _ and expand explanation
%                   'S11_dB' 
%                   'S11_complex' - ToDo?
%                   'Gen' - Ignored (default)
% TODO_DWW: fix comments up

if ( strcmp(R.t, goalResType)  )  % Same type, just return
    Ri = R;
elseif ( strncmp(R.t,'S',1) && strncmp(goalResType,R.t,3) )  % S-parameters, e.g. S11
    % Assume that the response is always complex for the response because inner workings handle it like this.
    s = regexp(goalResType, '\_', 'split');
    assert(length(s) == 2);
    
    switch s{2}
        case 'complex'
            Ri = R;
        case 'dB'
            Ri = R;
            Ri.r = dB20(R.r);
        otherwise
            % TODO_DWW: Flesh this out
            error(['Unknows combination found. R.t = ', R.t, ', goalResType = ', goalResType])
    end
% else if (  )        % General
% TODO_DWW: Flesh this out
else
    error(['Unknows combination found. R.t = ', R.t, ', goalResType = ', goalResType])
end

end % convertResponse function