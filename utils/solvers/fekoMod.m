function R = fekoMod(M, xi, xp, Rtype, f)
% TODO_DWW: Frequency is not changed from the outside using the model like other routines.
% TODO_DWW: Add functunality for coarse model xp

assert(isempty(xp), 'Implicit parameters has not been implemented for this routine.');
assert(isempty(f), 'Model independent frequency for coarse models has not been implemented yet.')
% Nq = length(xp);

Nn = length(xi);
Nr = length(Rtype);
R = cell(1,Nr);

% Run FEKO - cannot run with path, so change the directory
curDir = pwd;
cd(M.path);
% Create temp file for run so that the original isn't overwritten and doesn't show up in version control.
copyfile([M.name,'.cfx'], ['temp_',M.name,'.cfx'])
M.name = ['temp_',M.name];

% Build parameter string
parStr = [];
for nn = 1:Nn
    parStr = [parStr,' -#',M.params{nn},'=',num2str(xi(nn))];
end
% Remesh the structure with the new parameters
FEKOmesh = ['cadfeko_batch ',[M.path,M.name,'.cfx'],parStr];
[statusMesh, cmdoutMesh] = system(FEKOmesh);
% % Run FEKO - cannot run with path, so change the directory
% curDir = pwd;
% cd(M.path);
% % Create temp file for run so that the original isn't overwritten and doesn't show up in version control.
% copyfile([M.name,'.cfx'], ['temp_',M.name,'.cfx'])
% M.name = ['temp_',M.name];
FEKOrun = ['runfeko ', [M.name,'.cfx']];
[statusRun, cmdoutRun] = system(FEKOrun);
cd(curDir);
% Generate output
for rr = 1:Nr
    if strncmp(Rtype{rr},'S',1)
        portB = Rtype{rr}(2:find(Rtype{rr}==',')-1);
        portB_num = str2double(portB);
        portA = Rtype{rr}(find(Rtype{rr}==',')+1:end);
        portA_num = str2double(portA);
        % Read the Sb_a touchstone file - must be exported by the FEKO
        % file with the correct name - Name_S1_1.s*p!
        fileName = [M.name,'_S',portB,'_',portA,'.s*p'];
        wildPos = find(fileName=='*');
        fileDir = dir([M.path,fileName]);
        portCount = str2double(fileDir.name(wildPos));
        file = [fileDir.folder, '\', fileDir.name];
        [Spar,freq] = touchread(file,portCount);
        Sba = reshape(Spar(portB_num,portA_num,:),length(freq),1);
        R{rr}.f = freq;
    else
        error(['Unrecognised Rtype{rr} for FEKO run. Rtype{',rr,'} = ', Rtype{rr}])
    end
    R{rr}.r = Sba;
    R{rr}.t = Rtype{rr};
end

save(['fekoLog', datestr(datetime('now'), '_yyyymmdd_HHMM'), '.mat'], 'statusMesh', ...
        'cmdoutMesh', 'statusRun', 'cmdoutRun');

fid=fopen(['fekoLog', datestr(datetime('now'), '_yyyymmdd_HHMM'), '.txt'],'w');
% fprintf(fid, [ cmdoutMesh '\n\n' statusRun '\n\n' cmdoutRun]);
fprintf(fid, '%s \n\n', [cmdoutMesh]');
fprintf(fid, '%s \n\n', [statusRun]');
fprintf(fid, '%s \n\n', [cmdoutRun]');
fclose(fid);

end % fekoMod function
