function plotModels(plotModelFlag, plotGoalFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts)

% Plot the goals and the responses for the initial fine, coarse, 
% optimised surrogate and aligned surrogate.

% Inputs:
%   plotModelFlag:  If false then no plots are made.
%   plotGoalFlag:   If false then the goal values are not added to the plot.
%   itNum:      The current iteration number.
%   Rci:        Coarse model response.
%   Rfi:        Fine model response.
%   Rsi:        The surrogate response.
%   Rsai:       The aligned surrogate response.
%   OPTopts:    Optimization options, used for determining what the type of response is and goal values.

if plotModelFlag && valuesAreValid()
    figure()
    freq = [];
    if isfield(Rfi{1}{1},'f')
        freq = Rfi{1}{1}.f;
    else
        freq = 1:length(Rfi{1}{1}.r);
    end

    % Only plot one figure per response type, and put all the goals on the
    % same figure
    uniqueResultTypes = unique(OPTopts.goalResType);
    Nrt = length(uniqueResultTypes);
    responsePlotted = zeros(1,Nrt);  % Check if the current response has been plotted already
    
    Ng = length(OPTopts.goalType);
    for gg = 1:Ng
        goalType = OPTopts.goalType{gg};
        splits = regexp(OPTopts.goalResType{gg}, '\_', 'split');
%         assert(length(splits) == 2, ['Expecting two parts to the goalResType, found: ', cell2mat(splits), ' from: ', OPTopts.goalResType]);
        
        if (length(splits) == 2) && (strcmp(splits{2},'complex'))
            [~,plotRow] = ismember(OPTopts.goalResType(gg),uniqueResultTypes);
            
            subplot(Nrt, 2, plotRow*2-1)
            if ~(responsePlotted(plotRow))
                plotResponse([splits{1}, '_real']);
            end
            plotGoal(plotGoalFlag, OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, real(OPTopts.goalVal{gg}), goalType, freq)
            legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', ['imag',goalType])
            legend('Location','best')
            subplot(Nrt, 2, plotRow*2)
            if ~(responsePlotted(plotRow))
                plotResponse([splits{1}, '_imag']);
                responsePlotted(plotRow) = 1;
            end
            plotGoal(plotGoalFlag, OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, imag(OPTopts.goalVal{gg}), goalType, freq)
            legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', ['imag',goalType])
            legend('Location','best')
        else
            [~,plotRow] = ismember(OPTopts.goalResType(gg),uniqueResultTypes);
            subplot(Nrt, 2, [plotRow*2-1, plotRow*2])
            if ~(responsePlotted(plotRow))
                plotResponse(OPTopts.goalResType{gg});
                responsePlotted(plotRow) = 1;
            end
            if isfield(OPTopts, 'goalVal')
                plotGoal(plotGoalFlag, OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, OPTopts.goalVal{gg}, goalType, freq)
            end
            legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', goalType)
            legend('Location','best')
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
    [freqPlot, freqUnit] = freqScale(freq);
    
    plot(freqPlot, RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
    plot(freqPlot, RciMatch.r,'r','LineWidth',1.5), grid on, hold on
    plot(freqPlot, RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
    plot(freqPlot, RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
    xlabel(['Frequency (',freqUnit,')'])
else
    plot(RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
    plot(RciMatch.r,'r','LineWidth',1.5), grid on, hold on
    plot(RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
    plot(RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
    xlabel('Index')
end
ylabel(replace(goalResType,'_',' '))    % Remove underscores in the labels
title(['Iteration: ',num2str(itNum), ', goal of result type : ', replace(goalResType,'_',' ')])
end % plotForRealValuedResponse function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotGoal(plotGoalFlag, goalStart, goalStop, goalValue, goalType, freq)

if plotGoalFlag
    % Scale the frequencies for plotting
    [freqPlot, ~, scaleFact] = freqScale(freq);
    goalStartPlot = goalStart./scaleFact;
    goalStopPlot = goalStop./scaleFact;
    
    stdGoalPlot = 0;
    colour = 'k';
    switch goalType
    case 'lt'
        colour = 'c';
        stdGoalPlot = 1;
    case 'gt'
        colour = 'm';
        stdGoalPlot = 1;
    case 'eq'
        colour = 'g';
        stdGoalPlot = 1;
    case 'minimax'
        colour = 'y';
        stdGoalPlot = 1;
    case 'bw'
        colour = 'b';
    case 'nPeaks'
        % TODO: Add plots for nPeaks
    case 'peaksVal'    
        colour = 'y';
        stdGoalPlot = 1;
    otherwise
        error(['Unknows goalType found. OPTopts.goalType = ', goalType])
    end

    if ( isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') && stdGoalPlot)
        iValidFreq = freq > goalStart & freq < goalStop;
        goalFreqPoints = freqPlot(iValidFreq);
    %     goalFreqPoints = freq(freq > goalStart & freq < goalStop);
        if length(goalValue) == 1  
            plot([goalStartPlot, goalStopPlot], (goalValue)*ones(1,2), colour, 'LineWidth',3)
            plot(goalFreqPoints, (goalValue)*ones(1,length(goalFreqPoints)),['o',colour],'LineWidth',3)
        elseif length(goalValue) == length(freq)   
            plot(freqPlot, goalValue,colour,'LineWidth',3)
            plot(goalFreqPoints, goalValue(iValidFreq),['o',colour],'LineWidth',3)
        else
            error('Goal plotting error - goalVal should be either length 1 or length(freq)');
        end
    end
end % plotGoal flag

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

