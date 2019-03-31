function Z = S2Z(S,ZL)

% function Z = S2Z(S,ZL)
% converts the NxN S matrix to the NxN Z matrix, where Zl is the scalar 
% (or vector of length N) characteristic impedances at each port

% Dirk de Villiers
% Create: 2017-01-31
% Last edit: 2017-01-31

[N1,N2] = size(S);
if N1 ~= N2
    error('S must be square');
end
NZL = length(ZL);
if NZL < 1 && NZL ~= N1
    error('ZL should be scalar or of length N - the size of S');
end

if NZL == 1
    sz = eye(N1).*sqrt(ZL);
else
    sz = diag(sqrt(ZL));
end
% Z = sz*(eye(N1) + S)*inv(eye(N1) - S)*sz;
Z = sz*(eye(N1) + S)/(eye(N1) - S)*sz;


