function plotModels(plotModelFlag, plotGoalFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts, plotOpts)

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
%   plotOpts:   Graph options:
%       - legendLocation: defaults to 'Best'
%       - xlim: defaults to min and max of freqeuency.
%       - ylim: defaults to internal plotting defaults.
%       - xtick: defaults to internal plotting defaults.
%       - ytick: defaults to internal plotting defaults.
%       - pbaspectVec: defaults to internal plotting defaults.

if plotModelFlag
    
    figure()
    % ----- Configure frequencies -----  
    freq = [];
    if isfield(Rfi{1}{1},'f')
        freq = Rfi{1}{1}.f;
        [freqPlot, freqUnit, scaleFact] = freqScale(freq);
    else
        freqPlot = 1:length(Rfi{1}{1}.r);
        freqUnit = 1;
    end

    % ----- Set defaults ----- 
    legendLocation = 'Best';
    xlim = [min(freqPlot), max(freqPlot)];
    ylim = [];
    xtick = [];
    ytick = [];
    pbaspectVec = [];
    if isfield(plotOpts, 'legendLocation'), legendLocation = plotOpts.legendLocation; end
    if isfield(plotOpts, 'xlim'), xlim = plotOpts.xlim; end
    if isfield(plotOpts, 'ylim'), ylim = plotOpts.ylim; end
    if isfield(plotOpts, 'xtick'), xtick = plotOpts.xtick; end
    if isfield(plotOpts, 'ytick'), ytick = plotOpts.ytick; end
    if isfield(plotOpts, 'pbaspectVec'), pbaspectVec = plotOpts.pbaspectVec; end
            
    % ----- Do plotting -----  
    if isfield(OPTopts, 'goalResType') && valuesAreValid()
        detailedPlotting()
    else
        basicPlotting()
    end

end % if flag and validation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function detailedPlotting()
% Only plot one figure per response type, and put all the goals on the
% same figure
markerstr = 'xso+*d^v><ph.xso+*d^v><ph.';
colourstr = 'krbgmcykrbgmckrbgmcykrbgmc';
fontsize = 13;

uniqueResultTypes = unique(OPTopts.goalResType);
Nrt = length(uniqueResultTypes);
responsePlotted = zeros(1,Nrt);  % Check if the current response has been plotted already



legendText = {'Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate'};
% keyboard
Ng = length(OPTopts.goalType);
% TODO_DWW: DDV: No longer plots each goal seperately.
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
        plotGoal(plotGoalFlag, OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, real(OPTopts.goalVal{gg}), goalType, freqPlot, freqUnit, scaleFact)
        legendText{end+1} = ['imag',goalType];
        % legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', ['imag',goalType])
        % legend('Location',legendLocation)
        subplot(Nrt, 2, plotRow*2)
        if ~(responsePlotted(plotRow))
            plotResponse([splits{1}, '_imag']);
            responsePlotted(plotRow) = 1;
        end
        plotGoal(plotGoalFlag, OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, imag(OPTopts.goalVal{gg}), goalType, freqPlot, freqUnit, scaleFact)
        legendText{end+1} = ['imag',goalType];
        % legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate', )
        % legend('Location',legendLocation)
    else
        [~,plotRow] = ismember(OPTopts.goalResType(gg),uniqueResultTypes);
        subplot(Nrt, 2, [plotRow*2-1, plotRow*2])
        if ~(responsePlotted(plotRow))
            plotResponse(OPTopts.goalResType{gg});
            responsePlotted(plotRow) = 1;
        end
        if isfield(OPTopts, 'goalVal')
            plotGoal(plotGoalFlag, OPTopts.goalStart{gg}, OPTopts.goalStop{gg}, OPTopts.goalVal{gg}, goalType, freqPlot, freqUnit, scaleFact)
        end
        legendText{end+1} = goalType;
    end
end % for Ng

legend(legendText)
legend('Location',legendLocation)
set(gca, ...
'FontSize', fontsize);
if ~isempty(xlim), set(gca, 'xlim', xlim); end
if ~isempty(ylim), set(gca, 'ylim', ylim); end
if ~isempty(xtick), set(gca, 'XTick', xtick); end
if ~isempty(ytick), set(gca, 'YTick', ytick); end
if ~isempty(pbaspectVec), set(gca, 'PlotBoxAspectRatio', pbaspectVec); end

end % detailedPlotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basicPlotting()
markerstr = 'xso+*d^v><ph.xso+*d^v><ph.';
colourstr = 'krbgmcykrbgmckrbgmcykrbgmc';
fontsize = 13;

justMarkers = 1;

if ~isempty(Rfi) && isfield(Rfi{itNum}{1}, 't')
    RfiMatch = findResponseFor(Rfi{itNum}, Rfi{itNum}{1}.t);
    if isfield(RfiMatch,'f')
        freq = RfiMatch.f;
        [freqPlot, freqUnit] = freqScale(freq);
        if justMarkers 
            plot(freqPlot, RfiMatch.r, strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(freqPlot, RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
       
        xlabel(['Frequency (',freqUnit,')'])
    else
        if justMarkers 
            plot(RfiMatch.r, strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
        xlabel('Index')
    end
end

if ~isempty(Rci) && isfield(Rci{itNum}{1}, 't')
    RciMatch = findResponseFor(Rci{itNum}, Rci{itNum}{1}.t);
    if isfield(RciMatch,'f')
        freq = RciMatch.f;
        [freqPlot, freqUnit] = freqScale(freq);
        if justMarkers 
            plot(freqPlot, RciMatch.r, strcat(markerstr(2),colourstr(2)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(freqPlot, RciMatch.r,'r','LineWidth',1.5), grid on, hold on
        xlabel(['Frequency (',freqUnit,')'])
    else
        if justMarkers 
            plot(RciMatch.r,strcat(markerstr(2),colourstr(2)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(RciMatch.r,'r','LineWidth',1.5), grid on, hold on
        xlabel('Index')
    end
end

if ~isempty(Rsi) && isfield(Rsi{itNum}{1}, 't')
    RsiMatch = findResponseFor(Rsi{itNum}, Rsi{itNum}{1}.t);
    if isfield(RsiMatch,'f')
        freq = RsiMatch.f;
        [freqPlot, freqUnit] = freqScale(freq);
        if justMarkers 
            plot(freqPlot, RsiMatch.r, strcat(markerstr(3),colourstr(3)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(freqPlot, RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
        xlabel(['Frequency (',freqUnit,')'])
    else
        if justMarkers 
            plot(RsiMatch.r,strcat(markerstr(3),colourstr(3)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
        xlabel('Index')
    end
end

if ~isempty(Rsai) && isfield(Rsai{itNum}{1}, 't')
    RsaiMatch = findResponseFor(Rsai{itNum}, Rsai{itNum}{1}.t);
    if isfield(RsaiMatch,'f')
        freq = RsaiMatch.f;
        [freqPlot, freqUnit] = freqScale(freq);
        if justMarkers 
            plot(freqPlot, RsaiMatch.r, strcat(markerstr(4),colourstr(4)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(freqPlot, RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
        xlabel(['Frequency (',freqUnit,')'])
    else
        if justMarkers 
            plot(RsaiMatch.r,strcat(markerstr(4),colourstr(4)),'LineWidth',2,'MarkerSize',10), grid on, hold on
        else
        end
        plot(RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
        xlabel('Index')
    end
end

set(gca, ...
    'FontSize', fontsize);

% Assuming all responses are based off the same freq
if isfield(RciMatch,'f')
    freq = RciMatch.f;
    [freqPlot, freqUnit] = freqScale(freq);
    xlabel(['Frequency (',freqUnit,')'])
else
    xlabel('Index')
end

% ylabel(replace(goalResType,'_',' '))    % Remove underscores in the labels
% title(['Iteration: ',num2str(itNum), ', goal of result type : ', replace(goalResType,'_',' ')])
end % basicPlotting



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotResponse(goalResType)
% Find and convert responses to the particular goal type for plotting
RfiMatch = findResponseFor(Rfi{itNum}, goalResType);
RciMatch = findResponseFor(Rci{itNum}, goalResType);
RsiMatch = findResponseFor(Rsi{itNum}, goalResType);
RsaiMatch = findResponseFor(Rsai{itNum}, goalResType);

% Assuming all responses are based off the same freq
plot(freqPlot, RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
plot(freqPlot, RciMatch.r,'r','LineWidth',1.5), grid on, hold on
plot(freqPlot, RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
plot(freqPlot, RsaiMatch.r,'g--','LineWidth',2), grid on, hold on

if isfield(RciMatch,'f')
    % freq = RciMatch.f;
    % [freqPlot, freqUnit] = freqScale(freq);
    xlabel(['Frequency (',freqUnit,')'])
else
    % plot(RfiMatch.r,'k','LineWidth',1.5), grid on, hold on
    % plot(RciMatch.r,'r','LineWidth',1.5), grid on, hold on
    % plot(RsiMatch.r,'b:','LineWidth',1.5), grid on, hold on
    % plot(RsaiMatch.r,'g--','LineWidth',2), grid on, hold on
    xlabel('Index')
end
ylabel(replace(goalResType,'_',' '))    % Remove underscores in the labels
title(['Iteration: ',num2str(itNum), ', goal of result type : ', replace(goalResType,'_',' ')])

end % plotForRealValuedResponse function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotGoal(plotGoalFlag, goalStart, goalStop, goalValue, goalType, freq, freqUnit, scaleFact)

if plotGoalFlag
    % Scale the frequencies for plotting
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
case 'eqPhaseTune'
    colour = 'g';
    stdGoalPlot = 0;
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
    error(['Unknown goalType found. OPTopts.goalType = ', goalType])
end

    if ( isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') && stdGoalPlot)
        iValidFreq = freqPlot > goalStart & freqPlot < goalStop;
        goalFreqPoints = freqPlot(iValidFreq);
    %     goalFreqPoints = freq(freq > goalStart & freq < goalStop);
        if length(goalValue) == 1  
            plot([goalStartPlot, goalStopPlot], (goalValue)*ones(1,2), colour, 'LineWidth',3)
            plot(goalFreqPoints, (goalValue)*ones(1,length(goalFreqPoints)),['o',colour],'LineWidth',3)
        elseif length(goalValue) == length(freqPlot)   
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

