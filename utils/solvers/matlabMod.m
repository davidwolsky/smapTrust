function R = matlabMod(M, xi, xp, Rtype, f)

Nr = length(Rtype);
R = cell(1,Nr);

Ni = length(M.params);  % This is interpreted as the number of inputs to the function
inType = [];
for ii = 1:Ni
    inType = [inType, M.params{ii}];   % Initialise
end
switch inType
    case 'xf'
        Ri = M.name(xi, M.freq);
    case 'xxpf'
        Ri = M.name(xi, xp, M.freq);
    case 'xxp'
        Ri = M.name(xi, xp);
    otherwise
        Ri = M.name(xi);
end
% Distribute the responses
% For MATLAB case the model should the return the specified responses
% columnwise...
for rr = 1:Nr
    R{rr}.r = Ri(:,rr);
    R{rr}.t = Rtype{rr};	
    if ( exist('f','var') && ~isempty(f) )
        R{rr}.f = M.freq;
    end
end

end % matlabMod function