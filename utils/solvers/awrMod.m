function R = awrMod(M, xi, xp, Rtype, f)

% keyboard
% assert(~exist('xp','var'), 'This has only been implemented for coarse model evaluations...')
assert(logical(exist('xp','var')), 'This has only been implemented for coarse model evaluations...')
assert(isempty(f), 'Model independent frequency for coarse models has not been implemented yet.')
% TODO_DWW: Add functunality for fine model?

Nn = length(xi);
Nr = length(Rtype);
Nq = length(xp);
% TODO_DWW: If frequency is not specified then just take what is in the model already... like CST implementation.
Nm = length(M.freq);
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
    proj.Frequencies.Clear();
    for mm = 1:Nm
        awr.Project.Frequencies.Add(M.freq(mm));
    end
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
        [fin,Sxxin] = deal(zeros(nRead,1));
        fin = measurement_Sxx_Mag.XValues;

        for nn = 1:nRead
            amp = measurement_Sxx_Mag.YValue(nn,1);
            phase = measurement_Sxx_Ang.YValue(nn,1);
            Sxxin(nn) = amp.*exp(1i*phase);
        end

        if isfield(M,'freq')
            Nm = length(M.freq);
            Rreal = reshape(interp1(fin,real(Sxxin),M.freq,'spline'),Nm,1);
            Rimag = reshape(interp1(fin,imag(Sxxin),M.freq,'spline'),Nm,1);
            R{rr}.r = Rreal + 1i*Rimag;
            R{rr}.f = M.freq;
        else
            Nm = nRead;
            R{rr}.r = Sxxin;
            R{rr}.f = fin;
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