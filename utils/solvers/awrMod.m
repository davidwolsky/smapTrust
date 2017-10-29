function R = awrMod(M, xi, xp, Rtype, f)

% assert(~exist('xp','var'), 'This has only been implemented for coarse model evaluations...')
assert(logical(exist('xp','var')), 'This has only been implemented for coarse model evaluations...')
% assert(isempty(f), 'Model independent frequency for coarse models has not been implemented yet.')
% TODO_DWW: Add functunality for fine model?

if ~isempty(f)
    % Use the frequency passed in.
    fc = f;
elseif isfield(M,'freq') && ~isempty(M.freq)
    % Use frequency from the model.
    fc = M.freq;
else
    % No a run using frequency, relying on the model being set up.
    warning('AWR requires a frequency to run, hoping that the model is set up corrctly.')
    fc = [];
end

% assert(~isempty(fc), 'A frequeny is required one way or another.')


Nn = length(xi);
Nr = length(Rtype);
Nq = length(xp);
Nm = length(fc);
R = cell(1,Nr);

% Ensure that NI AWR is open and that all projects have been closed.
awr = actxserver('AWR.MWOffice');
% pause(0.001);
if ( awr.ProjectOpen() )
    proj = awr.Project;
else
    [M.path,M.name,'.emp']
    awr.invoke('Open', [M.path,M.name,'.emp']);
    proj = awr.Project;
end

if (proj.Frequencies.Count == 0)
    assert(~isempty(fc), 'A frequeny is required one way or another.')
    for mm = 1:Nm
        awr.Project.Frequencies.Add(fc(mm));
    end
elseif (~isempty(fc))

    % TODO_DWW: CRC_DDV: Either of these really should work, any ideas?
    % proj.Frequencies.Item(1).Value
    % proj.Frequencies.invoke('Item', 1)

    freqItem = [];
    mm = 1;
    % This is an attempt at an optimisation to not have to clear the frequency 
    % collection every time.
    % - 'Double' freqItem returned means that the frequency already exists.
    % - 'Interface.AWR_Design_Environment_12.IFrequency' means that a new frequency was added
    while isequal(class(freqItem),'double') && mm <= Nm
        freqItem = awr.Project.Frequencies.Add(fc(mm));
        mm = mm +1;
    end

    % If a new frequency was added then the integrity of the list is lost. 
    % We may be adding more frequencies to an old set. 
    if isequal(class(freqItem),'Interface.AWR_Design_Environment_12.IFrequency')
        proj.Frequencies.Clear();
        for mm = 1:Nm
            awr.Project.Frequencies.Add(fc(mm));
        end
    elseif isequal(class(freqItem),'double')
        % Do nothing.
    else
        error('An unkown frequency class type was returned.')
    end

    assert(proj.Frequencies.Count == Nm, 'There must the the same number of frequencies in AWR and our array.')
else
    % Do nothing, running with the model as is.    
end

eqns = proj.GlobalDefinitionDocuments.Item(1).Equations;

for nn = 1:Nn
    if ~eqns.Exists(M.params{nn}) 
        error(strcat('AWR parameter not found: ', M.params{nn}))
    end
    parStr = [M.params{nn},'=',num2str(xi(nn))];
    eqns.Item(M.params{nn}).Expression = parStr;
end
% Also include possible implicit parameters
for qq = 1:Nq
    if ~eqns.Exists(M.Iparams{qq}) 
        error(['AWR parameter not found: ', M.Iparams{qq}])
    end
    parStr = [M.Iparams{qq},'=',num2str(xp(qq))];
    eqns.Item(M.Iparams{qq}).Expression = parStr;
end

% Generate output
for rr = 1:Nr
    if strncmp(Rtype{rr},'S',1)
        assert(~contains(Rtype{rr},'_'), 'S-parameters are always complex therefore specifying units here does not make sense. The goal type can have units though.')
        RTypeString = replace(Rtype{rr},',','_');
        portsString = [Rtype{rr}(2:end)];

        % Adding a graph and measurement 
        graphs = proj.Graphs;
        if graphs.Exists([RTypeString,'_Mag Graph'])
            graph = graphs.Item([RTypeString,'_Mag Graph']);
        else
            graph = graphs.Add([RTypeString,'_Mag Graph'],'mwGT_Rectangular');
        end
        measurement_Sxx_Mag = graph.Measurements.Add(M.name,['|S(', portsString, ')|']);

        if graphs.Exists([RTypeString,'_Ang Graph'])
            graph = graphs.Item([RTypeString,'_Ang Graph']);
        else
            graph = graphs.Add([RTypeString,'_Ang Graph'],'mwGT_Rectangular');
        end
        measurement_Sxx_Ang = graph.Measurements.Add(M.name,['Ang(S(', portsString, '))']);

        proj.Simulator.Analyze;
        
        switch proj.ErrorState
            case 'mwET_Fatal'
                % R{rr}.hasError = true;
                error(['AWR proj.ErrorState gave mwET_Fatal.'])
            case 'mwET_Error'
                % R{rr}.hasError = true;
                error(['AWR proj.ErrorState gave mwET_Error.'])
            case 'mwET_Warning'
                % R{rr}.hasError = true;
                error(['AWR proj.ErrorState gave mwET_Warning.'])
            case 'mwET_None'
                % R{rr}.hasError = false;
                % Do nothing, continue
            otherwise
                % R{rr}.hasError = true;
                error(['AWR proj.ErrorState unknown result from simulation.'])
        end

        nRead = measurement_Sxx_Mag.XPointCount;
        [xin,Sxxin] = deal(zeros(nRead,1));
        xin = measurement_Sxx_Mag.XValues;

        for nn = 1:nRead
            amp = measurement_Sxx_Mag.YValue(nn,1);
            phase = measurement_Sxx_Ang.YValue(nn,1);
            Sxxin(nn) = amp.*exp(1i*phase);
        end

        if (~isempty(fc))
            % If a frequency has been specified use it...
            Nm = length(M.freq);
            Rreal = reshape(real(Sxxin), Nm, 1);
            Rimag = reshape(imag(Sxxin), Nm, 1);
            R{rr}.r = Rreal + 1i*Rimag;
            R{rr}.f = fc;

            assert(isequal(reshape(xin,Nm,1),fc), 'The xaxis and the freq requested should match.')
            % Would use interpolation on value is the frequency were to not change.
            % Rreal = reshape(interp1(xin,real(Sxxin),fc,'spline'),Nm,1);
            % Rimag = reshape(interp1(xin,imag(Sxxin),fc,'spline'),Nm,1);
        else
            % Nm = nRead;
            R{rr}.r = Sxxin;
            R{rr}.f = xin;
        end
        R{rr}.t = Rtype{rr};
    else
        error(['Unrecognised Rtype request ', Rtype{rr}]);
    end
end % for Nr
% awr.Project.Close(false)
% awr.Quit()
release(awr)

end % awrMod function