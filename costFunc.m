function cost = costFunc(R,GOALS)

% Function to calculate the cost of a given response R as specified by the
% goals in GOALS
% R is a cell array of structures containing the response in R.r, the type R.t, and the
% (optional) domain (typically frequency) in R.f.
% R can also be a structure if only one type of response is considered.
% If no response type or domain is specified, R may also be a vector.
% GOALS can contain (typically a subset of OPTopts used in the main function):
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'S11_dB' 
%               'S11_complex' - ToDo?
%               'Gen' - Ignored (default)
%   goalType:   Cell array of goal types {1,Ng}
%               Valid types:
%                   'lt' (Less than)
%                   'gt' (Greater than)
%                   'eq' (equal to)
%                   'minimax'
%                   'bw' (Like 'lt' but also maximizes bandwidth - requires goalCent value)
%   goalVal:    Cell array of goal values {1,Ng} - same order as goalType
%   goalWeight: Vector of goal weights [1,Ng] (default equal weights)
%   goalStart:  Cell array of start of valid goal domain {1,Ng} (optional)
%   goalStop:   Cell array of stop of valid goal domain {1,Ng} (optional)
%   goalCent:   Cell array of centre point of goal domain {1,Ng} (used by the 'bw' goalType) (optional)
%   errNorm:    Cell array of type of error norms to use for optimization {1,Ng}
%               Valid types: integers or inf
%                   
% Date created: 2015-06-26
% Dirk de Villiers 
% Last Modified: 2015-06-26
% Updates:
% 2015-06-26: Write function shell and basic functionality from SMmain.m migration


% Make R a structure if only a vector is passed
if ~isstruct(R) && ~iscell(R), R.r = R; end;
% Make R a cell array if only one structure is passed.
if length(R) == 1 && ~iscell(R), R = {R}; end

Nr = length(R);
Ng = length(GOALS.goalType);
[cSum,wSum] = deal(0);

for gg = 1:Ng
    G.goalResType = 'Gen';
    if isfield(GOALS,'goalResType'),G.goalResType = GOALS.goalResType{gg}; end
    G.goalType = GOALS.goalType{gg};
    if isfield(GOALS,'goalVal'),G.goalVal = GOALS.goalVal{gg}; end
    if isfield(GOALS,'goalStart'), G.goalStart = GOALS.goalStart{gg}; end
    if isfield(GOALS,'goalStop'), G.goalStop = GOALS.goalStop{gg}; end
    if isfield(GOALS,'goalCent'), G.goalCent = GOALS.goalCent{gg}; end
    G.errNorm = 'L1';
    if isfield(GOALS,'errNorm'), G.errNorm = GOALS.errNorm{gg}; end
    G.goalWeight = 1;
    if isfield(GOALS,'goalWeight'), G.goalWeight = GOALS.goalWeight{gg}; end
    wSum = wSum + G.goalWeight;
    
    foundMathingType = false;
    tt = 1;
    % TODO_DWW: Move this finding loop to a function and share with plotModels
    while tt <= Nr
        % Special case for complex and dB S11 goals...
        % TODO_DWW: Clean up -> I dont think we are going to handle this case any more.
        % if isfield(R{tt},'t') && isfield(GOALS,'goalResType') && strcmp(R{tt}.t,'S11') && strcmp(G.goalResType,'S11_dB')
        if isfield(R{tt},'t') && isfield(GOALS,'goalResType') && strncmp(GOALS.goalResType{gg},R{tt}.t,3)
            Ri = convertResponse(R{tt}, GOALS.goalResType{gg});
            % Ri.r = dB20(R{tt}.r);
            foundMathingType = true;
            break;
        end
        % else if ( strncmp(R{tt}.t,'S11',3) && strncmp(GOALS.goalResType{gg},'S11',3) ) 
        % if strcmp(R{tt}.t,GOALS.goalResType{gg})
            % Ri = R{tt};
            % foundMathingType = true;
            % break;
        % end
        tt = tt + 1;
    end
    assert(foundMathingType, ['No matching result type was found for the specified goalResType.', GOALS.goalResType{gg}, '.  R{:}.t = ', R{:}.t])
    % assert(isfield(GOALS,'goalResType') && ~strcmp(G.goalResType,'S11_dB'), 'G.goalResType = S11_dB is not supported at this time.')

    % Sort out the centre, start and stop positions - use index if no frequency is
    % specified in the response structure
    Nm = length(Ri.r);
    iStart = 1;
    iStop = Nm;
    iCent = round(Nm/2);
    if ~isfield(Ri,'f')
        if isfield(G,'goalStart')
            iStart = G.goalStart;
        end
        if isfield(G,'goalStop')
            iStop = G.goalStop;
        end
        if isfield(G,'goalCent')
            iCent = G.goalCent;
        end
    else
        if isfield(G,'goalStart')
            iStart = find(Ri.f >= G.goalStart,1);
        end
        if isfield(G,'goalStop')
            iStop = find(Ri.f <= G.goalStop,1,'last');
        end
        if isfield(G,'goalCent')
            iCent = find(Ri.f <= G.goalCent,1,'last');
        end
    end

    % Get the goal type and calculate the cost
    Rvalid = Ri.r(iStart:iStop);
    switch G.goalType
        case 'lt'
            y = Rvalid - G.goalVal;
            y(y < 0) = 0;
            c0 = norm(y,G.errNorm);
            
        case 'gt'
            y = Rvalid - G.goalVal;
            y(y > 0) = 0;
            c0 = norm(y,G.errNorm);
            
        case 'eq'
            if length(G.goalVal) == length(Ri.r)
                y = Rvalid - G.goalVal(iStart:iStop);
            else
                y = Rvalid - G.goalVal;
            end
            c0 = norm(y,G.errNorm);
            
        case 'minimax'
            c0 = max(Rvalid);
            
        case 'bw'
            y = Rvalid - G.goalVal;
            iValid = find(y < 0);
            if isempty(iValid)
                c0 = Nm + min(y);
            else
                % Get the lowest in band index
                i1 = iValid(1);
                % Get the highest in band index
                i2 = iValid(end);
                % Get number of indexes to estimate a sensible penalty factor
                Ni = i2-i1;
                % New valid region response
                yi = y(i1:i2);
                yi(yi < 0) = 0;
                if max(yi) == 0
                    beta = 0;
                else
                    beta = 10*Ni./max(yi).^2;
                end
%                 c0 = -min(i2-iCent,iCent-i1) +  beta*Lnorm(yi,'L2');
                  c0 = -min(i2-iCent,iCent-i1) +  beta*norm(yi,2);
            end
        otherwise
            error('Unknown goalType');
    end
    
    cSum = cSum + G.goalWeight*c0;
end
cost = cSum/wSum;

end