function plotModels(plotFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts)

% Plot the goals and the responses for the initial fine, coarse, 
% optimised surrogate and aligned surrogate.

% Inputs:
%   plotFlag:   If false then no plots are made.
%   itNum:      The current iteration number.
%   Rci:        Coarse model response.
%   Rfi:        Fine model response.
%   Rsi:        The surrogate response.
%   Rsai:       The alligned surrogate response.
%   OPTopts:    Operational options, used for determining what the type of response is and goal values.

if plotFlag && valuesAreValid()
    figure()
    freq = [];
    if isfield(Rfi{1}{1},'f')
        freq = Rfi{1}{1}.f;
    else
        freq = 1:length(Rfi{1}{1}.r);
    end

    Ng = length(OPTopts.goalType);
    for gg = 1:Ng
        goalType = OPTopts.goalType{gg};
        splits = regexp(OPTopts.goalResType{gg}, '\_', 'split');
%         assert(length(splits) == 2, ['Expecting two parts to the goalResType, found: ', cell2mat(splits), ' from: ', OPTopts.goalResType]);
        
        if (length(splits) == 2) && (strcmp(splits{2},'complex'))
            subplot(Ng, 2, gg*2-1)
            plotResponse([splits{1}, '_real']);
            plotGoal(OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, real(OPTopts.goalVal{gg}), goalType, freq)
            subplot(Ng, 2, gg*2)
            % legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', ['imag',goalType])
            plotResponse([splits{1}, '_imag']);
            plotGoal(OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, imag(OPTopts.goalVal{gg}), goalType, freq)
            
            % legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', ['imag',goalType])
        else
            subplot(Ng, 2, [gg*2-1, gg*2])
            plotResponse(OPTopts.goalResType{gg});
            if isfield(OPTopts, 'goalVal')
                plotGoal(OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, OPTopts.goalVal{gg}, goalType, freq)
            end
            
            % legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', goalType)
        end
    end % for Ng
end % if validation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotResponse(goalResType)
% Find and convert responses to the particular goal type for plotting
RfiMatch = findResponseFor(Rfi{itNum}, goalResType);
RciMatch = findResponseFor(Rci{itNum}, goalResType);
RsiMatch = findResponseFor(Rsi{itNum}, goalResType);
RsaiMatch = findResponseFor(Rsai{itNum}, goalResType);

% Assuming all responses are based off the same freq
if isfield(RciMatch,'f')
    freq = RciMatch.f;
    plot(freq, RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
    plot(freq, RciMatch.r,'r','LineWidth',1.5), grid on, hold on
    plot(freq, RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
    plot(freq, RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
    xlabel('Frequency')
else
    plot(RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
    plot(RciMatch.r,'r','LineWidth',1.5), grid on, hold on
    plot(RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
    plot(RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
    xlabel('Index')
end
ylabel(goalResType)
title(['Iteration: ',num2str(itNum), ', goal of result type : ', goalResType])
end % plotForRealValuedResponse function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotGoal(goalStart, goalStop, goalValue, goalType, freq)
switch goalType
case 'lt'
    colour = 'c';
case 'gt'
    colour = 'm';
case 'eq'
    colour = 'g';
case 'minimax'
    colour = 'y';
case 'bw'
    colour = 'b';
otherwise
    error(['Unknows goalType found. OPTopts.goalType = ', goalType])
end

if ( isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    plot([goalStart, goalStop], (goalValue)*ones(1,2), colour, 'LineWidth',3)
    
    goalFreqPoints = freq(freq > goalStart & freq < goalStop);
    plot(goalFreqPoints, (goalValue)*ones(1,length(goalFreqPoints)),['o',colour],'LineWidth',3)
end
end % plotGoal function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [areValid] = valuesAreValid()
areValid = true;
areValid = areValid & isfield(OPTopts,'goalResType');
areValid = areValid & (length(Rci) >= itNum);
areValid = areValid & (length(Rfi) >= itNum);
areValid = areValid & (length(Rsi) >= itNum);
areValid = areValid & (length(Rsai) >= itNum);
end % valuesAreValid function

end % plotModels function
