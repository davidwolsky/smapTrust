function Res = convertResponse(Ri, goalResType)

% Function to convert a response (R) from one type (R.t) to a specified goalResType. 
% Ri is a single instance of structures containing the response in Ri.r, the type Ri.t, and the
% (optional) domain (typically frequency) in R.f.
% Ri can also be a structure if only one type of response is considered.
% If no response type or domain is specified, Ri may also be a vector.
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

if ( strcmp(Ri.t, goalResType)  )  % Same type, just return
    Res = Ri;
elseif ( strncmp(goalResType,Ri.t,find(goalResType=='_')-1) )
    % Assume that the response is always complex for the response, because 
    % inner workings handle it like this.
    switch extractUnitString(goalResType)
        case 'complex'
            Res = Ri;
        case 'real'
            Res = Ri;
            Res.r = real(Ri.r);
        case 'imag'
            Res = Ri;
            Res.r = imag(Ri.r);
        case 'dB'
            Res = Ri;
            Res.r = dB20(Ri.r);
        case 'abs'
            Res = Ri;
            Res.r = abs(Ri.r);
        case 'angle'
            Res = Ri;
            Res.r = angle(Ri.r);
        case 'deg'
            Res = Ri;
            Res.r = deg(angle(Ri.r));
        otherwise
            error(['Unknows combination found. Ri.t = ', Ri.t, ', goalResType = ', goalResType])
    end
else
    error(['Unknows combination found. R.t = ', Ri.t, ', goalResType = ', goalResType])
end

end % convertResponse function