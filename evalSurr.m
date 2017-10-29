function [Rs] = evalSurr(x, S)

% function [Rs] = evalSurr(x,S)
% Evaluates the surrogate model, S, at the input positions in the vector x,
% and returns the responses in Rs
% Based on the 2006 MTT paper by Koziel on a Generalized SM Framework
%
% Returns:
% Rs returns Nm responces (like for m frequency points) for each input 
%      point as a matrix [Nm,Nc]
%
% Inputs:
% x is the matrix [Nn,Nc] of input parameters, each column corresponding to a
% different point in the parameter space.  Number of rows are the number
% of input parameters accepted by the models = Nn.
%
% The surrogate model is a structure containing:
% coarse:  Function handle pointing to a predefined function that evaluates
%          the coarse model operating as Rc = coarse(xc,xp,f) thus returning Rc 
%          (coarse model response) for each element in xc. Optional xp are
%          the implicit (pre-assigned) parameters, and optional f the
%          positions in the frequency domain where the model is evaluated.
%          If f is specified and there are no implicit parameters, xp must
%          be passed as an empty matrix.
%       Rc [Nm,Nc]
%       xc [Nn,Nc]
%       xp [Nq,1]
%       f [Nm,1]
% M:       Structure containing auxiliary parameters to pass to coarse.  See SMmain.m for details. 
% A:       Multiplicative OSM factor diag[Nm,Nm]
% B:       Multiplicative input SM factor [Nn,Nn]
% c:       Additative input SM term [Nn,1]
% xp:      Pre-assigned parameters [Nq,1]
% G:       Multiplicative ISM factor [Nq,Nn]
% d:       Additive zeroth order OSM term [Nm,1]
% E:       Additive first order OSM term [Nm,Nn]
% xi:      Position of last update point of the model [Nn,1]
% F:       Frequency space mapping parameters [2,1]
% f:       Coarse model frequency range [Nm,1]


% Date created: 2014-11-09
% Dirk de Villiers
% Last Modified: 2015-03-12
% Updates:
% 2014-11-09: Write function shell and basic functionality
% 2015-03-09: Start FSM - make shell
% 2015-03-10: Finish FSM 
%             Completely redefine the p-loop - simpler to understand now
%             Include M in S.  Used in the general main SM loop.
% 2015-03-12: Include limits on x and xp
% 2015-03-14: Remove (comment) limits on x and xp - should be handled externally...
% 2015-10-16: Improve updating of Rc in the Np loop - columns added
%             explicitly to be compatible with legacy case also

% Get vector sizes 
[Nn,Nc] = size(x);

% Set up defaults for the input SM 
if ~isfield(S, 'B')
    B = eye(Nn);
else
    B = S.B;
end
if ~isfield(S, 'c')
    c = zeros(Nn, Nc);
else
    c = repmat(S.c, 1, Nc);   % Repeat for each input point
end

% Check if any pre-assigned variables are provided and set up defaults for
% ISM
if isfield(S,'xp') && ~isempty(S.xp) 
    xp = S.xp;
    [Nq,dummy] = size(xp);
else
    xpc = [];
    Nq = 0;
end
if ~isfield(S,'G')
    G = zeros(Nq, Nn);
else
    G = S.G;
end

% Check for FSM and set up new frequency vector
F = [1;0];
if isfield(S, 'F'), F = S.F; end
if isfield(S, 'f')
    % A local frquency axis is used here. The result is returned without
    % the frequency so this function can indeed be used to observe the
    % results with a shift and scaling. The values are just slotted into 
    % array placeholders. 
    % 
    fc = F(1).*S.f + F(2);
else
    fc = [];
%     fc = S.f;
end

% Evaluate the coarse model to find the response vector size
Rc = [];
for cc = 1:Nc
    xc = B*x(:,cc) + c(:,cc);
    if isfield(S,'M')   % This will be the typical case - the else is left for legacy
        if Nq > 0
            xpc = G*x(:,cc) + xp;
        end
        RcStruct = S.coarse(S.M, xc, xpc, fc);
        Rc1 = RcStruct{1}.r;    % Assume this is returned from the SMmain.m function in the response structure format
    else    % Legacy case with no M (probably not going to be updated much...)
        if Nq==0 && ~isfield(S,'f')
            Rc1 = S.coarse(xc);
        elseif Nq==0 && isfield(S,'f')
            Rc1 = S.coarse(xc, [], fc);
        elseif ~isfield(S,'f')
            xpc = G*x(:,cc) + xp;
            Rc1 = S.coarse(xc, xpc);
        else
            xpc = G*x(:,cc) + xp;
            Rc1 = S.coarse(xc, xpc, fc);
        end
        Rc1 = Rc1;
%         Rc1 = reshape(Rc1,length(Rc1),1);
    end

    % The coarse models are appended so that output space mapping can be applied.
    % It is handled elsewhere to try and manage multiple reponse types for now... 
    Rc = [Rc, Rc1];
%     Rc(:,pp) = Rc1;
end

% Set up OSM defaults
[Nm,Nc] = size(Rc);
if ~isfield(S,'A')
    A = eye(Nm);
else
    A = S.A;
end
if ~isfield(S,'d')
    d = zeros(Nm,Nc);
else
    d = S.d;
    d = repmat(d,1,Nc); % repeat for all input points
end
if ~isfield(S,'E')
    E = zeros(Nm,Nn);
else
    E = S.E;
    if ~isfield(S,'xi')
        E = zeros(Nm,Nn);
        warning('E term ignored since no xi specified in S');
    else 
        
    end
end

if isfield(S,'xi')
    xi = repmat(S.xi,1,Nc); % Repeat for each input point
    xi = repmat(xi,Nm,1);   % And repeat for every response point
else    % This is just a dummy, since E will always be set to zero in this case...
    xi = x;
end

try
    Rs = A*Rc + d + E*(x-xi);
catch mE
    keyboard;
end


end % evalSurr function
