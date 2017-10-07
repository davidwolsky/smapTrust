function R = cstMod(M, xi, xp, Rtype, f)

% TODO_DWW: Add functunality for coarse model xp
assert(isempty(xp), 'Implicit parameters has not been implemented for this routine.')
% Nq = length(xp);
assert(isempty(f), 'Model independent frequency for coarse models has not been implemented yet.')

Nm = length(M.freq);

Nn = length(xi);
Nr = length(Rtype);
R = cell(1,Nr);

% Start CST activeX
cst = actxserver('CSTSTUDIO.Application');

% Get handle to the model
persistent mws
% Only open if the handle is empty - therefore on the first run
if isempty(mws)
    pause(0.1);
    mws = invoke(cst,'OpenFile',[M.path,M.name,'.cst']);
end

% Update parameters
for nn = 1:Nn
    invoke(mws,'StoreParameter',M.params{nn},xi(nn));
end
invoke(mws,'Rebuild');

% Get frequency unit
units = invoke(mws,'Units');
fUnit = invoke(units,'GetFrequencyUnit');
switch fUnit
    case 'KHz'
        fScale = 1e3;
    case 'MHz'
        fScale = 1e6;
    case 'GHz'
        fScale = 1e9;
    case 'THz'
        fScale = 1e12;
    case 'PHz'
        fScale = 1e15;
    otherwise
        fScale = 1;
end

% Run simulation

% TODO_DWW: Break up the CST flag to include F or T solver
solver = invoke(mws,'Solver');
% solver = invoke(mws,'FDSolver');      % This is the FD solver
invoke(solver,'Start');

% Generate output
for rr = 1:Nr
    if strncmp(Rtype{rr},'S',1)
        restree = invoke(mws,'Resulttree'); % Handle to the result tree
        SmnHandle = invoke(restree,'GetResultFromTreeItem',['1D Results\S-Parameters\',Rtype{rr}],'3D:RunID:0'); % Get the object handle
        nRead = invoke(SmnHandle,'GetN');   % How many frequency points?
        [fin,Smnre,Smnim] = deal(zeros(nRead,1));
        for nn = 1:nRead
            fin(nn) = invoke(SmnHandle,'GetX',nn-1);
            Smnre(nn) = invoke(SmnHandle,'GetYRe',nn-1);
            Smnim(nn) = invoke(SmnHandle,'GetYIm',nn-1);
        end
        finHz = fin.*fScale;    % Scaled to Hz
        Sin = Smnre + 1i.*Smnim;
        
        if isfield(M,'freq')
            Nm = length(M.freq);
            R{rr}.r = reshape(interp1(finHz,Sin,M.freq,'linear','extrap'),Nm,1);
            R{rr}.f = M.freq;
        else
            Nm = nRead;
            R{rr}.r = Sin;
            R{rr}.f = finHz;
        end
        R{rr}.t = Rtype{rr};
        release(SmnHandle);
    else
        error(['Only S-paramater type results supported in cstMod - I found: ', Rtype{rr},', which should start with an S if you want to stand any chance at success...']);
    end
end
invoke(mws,'Save');
% invoke(mws,'Quit');

% keyboard;

end % cstMod