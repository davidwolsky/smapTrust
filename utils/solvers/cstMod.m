function R = cstMod(M, xi, xp, Rtype, f)

% TODO_DWW: Add functunality for coarse model xp
assert(logical(~exist('xp','var')), 'Implicit parameters has not been implemented for this routine.')
% Nq = length(xp);
assert(isempty(f), 'Model independent frequency for coarse models has not been implemented yet.')

Nm = length(M.freq);

Nn = length(xi);
Nr = length(Rtype);
R = cell(1,Nr);

% Start CST activeX
cst = actxserver('CSTSTUDIO.Application');
% Get handle to the model
mws = invoke(cst,'OpenFile',[M.path,M.name,'.cst']);
% Update parameters
for nn = 1:Nn
    invoke(mws,'StoreParameter',M.params{nn},xi(nn));
end
invoke(mws,'Rebuild');
% Run simulation
solver = invoke(mws,'Solver');
invoke(solver,'Start');

% Generate output
for rr = 1:Nr
    % TODO_DWW: This needs to be updated
    if strcmp(Rtype{rr},'S1,1_dB')
        result = invoke(mws,'Result1D','d1(1)1(1)');    % Sb,a in dB
        % Get nr of frequency points in the plot
        nRead = invoke(result,'GetN');
        [fin,Sbain] = deal(zeros(nRead,1));
        for nn = 1:nRead
            fin(nn) = invoke(result,'GetX',nn-1);        % Typically in GHz
            Sbain(nn) = invoke(result,'GetY',nn-1);
        end
        if isfield(M,'freq')
            Nm = length(M.freq);
            R{rr}.r = reshape(interp1(fin,Sbain,M.freq,'spline'),Nm,1);
            R{rr}.f = M.freq;
        else
            Nm = nRead;
            R{rr}.r = Sbain;
            R{rr}.f = fin;
        end
        R{rr}.t = Rtype{rr};
        release(result);
    elseif strcmp(Rtype{rr},'S1,1_complex')
        resultA = invoke(mws,'Result1D','a1(1)1(1)');    % amplitude of Sb,a
        resultP = invoke(mws,'Result1D','p1(1)1(1)');    % phase of Sb,a
        % Get nr of frequency points in the plots
        nRead = invoke(resultA,'GetN');
        [fin,Sbain] = deal(zeros(nRead,1));
        for nn = 1:nRead
            fin(nn) = invoke(resultA,'GetX',nn-1);        % Typically in GHz
            amp = invoke(resultA,'GetY',nn-1);
            phase = rad(invoke(resultP,'GetY',nn-1));
            Sbain(nn) = amp.*exp(1i*phase);
        end
        if isfield(M,'freq')
            Nm = length(M.freq);
            Rreal = reshape(interp1(fin,real(Sbain),M.freq,'spline'),Nm,1);
            Rimag = reshape(interp1(fin,imag(Sbain),M.freq,'spline'),Nm,1);
            R{rr}.r = Rreal + 1i*Rimag;
            R{rr}.f = M.freq;
        else
            Nm = nRead;
            R{rr}.r = Sbain;
            R{rr}.f = fin;
        end
        R{rr}.t = Rtype{rr};
        release(resultA);
        release(resultP);
    end
end
invoke(mws,'Save');
invoke(mws,'Quit');


end % cstMod