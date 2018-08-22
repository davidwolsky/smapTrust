function R = matlabMod(M, xi, xp, Rtype, f)

Nr = length(Rtype);
R = cell(1,Nr);

if ~isempty(f)
    % Use the frequency passed in.
    fc = f;
elseif isfield(M,'freq') && ~isempty(M.freq)
    % Use frequency from the model.
    fc = M.freq;
else
    fc = [];
end

Ni = length(M.params);  % This is interpreted as the number of inputs to the function
inType = [];
for ii = 1:Ni
    inType = [inType, M.params{ii}];   % Initialise
end
switch inType
    case 'xf'
        Ri = M.name(xi, fc);
    case 'xxpf'
        Ri = M.name(xi, xp, fc);
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
    if ( exist('fc','var') && ~isempty(fc) )
        R{rr}.f = M.freq;
    end
end

end % matlabMod function