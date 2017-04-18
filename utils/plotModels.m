function plotModels(plotFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts)

% Plot the initial fine, coarse, optimised surrogate and aligned surrogate

% Inputs:
%   plotFlag:   If false then no plots are made.
%   itNum:      The current iteration number.
%   Rci:        Coarse model response.
%   Rfi:        Fine model response.
%   Rsi:        The surrogate response.
%   Rsai:       The alligned surrogate response.
%   OPTopts:    Operational options, used for determining what the type of response is and goal values.

if plotFlag && valuesAreValid()
    figure(itNum)
    
    % TODO_DWW: Take loop from costFunc to actually sort this nicely with goals associated with a response type
    %           plotted on the same graph and that the correct form of the R.t is plotted. dB and not complex...
	Nr = length(OPTopts.Rtype); % Number of responses requested
    if ( ~checkForComplexGoalResType() )
        columnCount = 1;
    else
        columnCount = 2;
    end

    % TODO_DWW: Add dB for goal type.

    for rr = 1:Nr
        subplot(Nr, columnCount, columnCount*(rr-1) + 1)
        plotForRealValuedResponse();
        if ( columnCount == 2 )
            subplot(Nr, columnCount, columnCount*(rr-1) + 2)
            plotForComplexResponse();
        end
    end % for Nr
end % if validation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotForRealValuedResponse()
if isfield(OPTopts,'goalResType') && strcmp(OPTopts.goalResType,'S11_dB')
    if isfield(Rci{itNum}{rr},'f')
        plot(Rci{itNum}{rr}.f,dB20(Rfi{itNum}{rr}.r),'k','LineWidth',1.5), grid on, hold on
        plot(Rci{itNum}{rr}.f,dB20(Rci{itNum}{rr}.r),'r','LineWidth',1.5), grid on, hold on
        plot(Rsi{itNum}{rr}.f,dB20(Rsi{itNum}{rr}.r),'b','LineWidth',1.5), grid on, hold on
        plot(Rsai{itNum}{rr}.f,dB20(Rsai{itNum}{rr}.r),'g--','LineWidth',2), grid on, hold on
        xlabel('Frequency')
    else
        plot(dB20(Rfi{itNum}{rr}.r),'k','LineWidth',1.5), grid on, hold on
        plot(dB20(Rci{itNum}{rr}.r),'r','LineWidth',1.5), grid on, hold on
        plot(dB20(Rsi{itNum}{rr}.r),'b','LineWidth',1.5), grid on, hold on
        plot(dB20(Rsai{itNum}{rr}.r),'g','LineWidth',1.5), grid on, hold on
        xlabel('Index')
    end
else
    if isfield(Rci{itNum}{rr},'f')
        plot(Rci{itNum}{rr}.f,real(Rfi{itNum}{rr}.r),'k','LineWidth',1.5), grid on, hold on
        plot(Rci{itNum}{rr}.f,real(Rci{itNum}{rr}.r),'r','LineWidth',1.5), grid on, hold on
        plot(Rsi{itNum}{rr}.f,real(Rsi{itNum}{rr}.r),'b','LineWidth',1.5), grid on, hold on
        plot(Rsai{itNum}{rr}.f,real(Rsai{itNum}{rr}.r),'g--','LineWidth',2), grid on, hold on
        xlabel('Frequency')
    else
        plot(real(Rfi{itNum}{rr}.r),'k','LineWidth',1.5), grid on, hold on
        plot(real(Rci{itNum}{rr}.r),'r','LineWidth',1.5), grid on, hold on
        plot(real(Rsi{itNum}{rr}.r),'b','LineWidth',1.5), grid on, hold on
        plot(real(Rsai{itNum}{rr}.r),'g','LineWidth',1.5), grid on, hold on
        xlabel('Index')
    end
end
plotRealGoals()
ylabel(OPTopts.Rtype{rr})
title(['Iteration ',num2str(itNum)])
legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate')
end % plotForRealValuedResponse function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotForComplexResponse()
if isfield(Rci{itNum}{rr},'f')
    plot(Rci{itNum}{rr}.f,imag(Rfi{itNum}{rr}.r),'k','LineWidth',1.5), grid on, hold on
    plot(Rci{itNum}{rr}.f,imag(Rci{itNum}{rr}.r),'r','LineWidth',1.5), grid on, hold on
    plot(Rsi{itNum}{rr}.f,imag(Rsi{itNum}{rr}.r),'b','LineWidth',1.5), grid on, hold on
    plot(Rsai{itNum}{rr}.f,imag(Rsai{itNum}{rr}.r),'g--','LineWidth',2), grid on, hold on
    xlabel('Frequency')
else
    plot(imag(Rfi{itNum}{rr}.r),'k','LineWidth',1.5), grid on, hold on
    plot(imag(Rci{itNum}{rr}.r),'r','LineWidth',1.5), grid on, hold on
    plot(imag(Rsi{itNum}{rr}.r),'b','LineWidth',1.5), grid on, hold on
    plot(imag(Rsai{itNum}{rr}.r),'g','LineWidth',1.5), grid on, hold on
    xlabel('Index')
end
plotImagGoals()
ylabel(OPTopts.Rtype{rr})
% TODO_DWW: update labels to reflect actual plot type, e.g. dB -> from goalResType
title(['Iteration ',num2str(itNum), ' - Imag part'])
legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate')
end % plotForComplexResponse function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRealGoals()
if ( isfield(OPTopts, 'goalVal') && isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        plot([OPTopts.goalStart{gg}, OPTopts.goalStop{gg}], real(OPTopts.goalVal{gg})*ones(1,2), 'm', 'LineWidth',3)
    end
end % if validation
end % plotRealGoals function
function plotImagGoals()
if ( isfield(OPTopts, 'goalVal') && isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        if ( ~isreal(OPTopts.goalVal{gg}) )
            plot([OPTopts.goalStart{gg}, OPTopts.goalStop{gg}], imag(OPTopts.goalVal{gg})*ones(1,2), 'm', 'LineWidth',3)
        end % imag part
    end
end % if validation
end % plotImagGoals function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [areValid] = valuesAreValid()
areValid = true;
areValid = areValid & (length(Rci) >= itNum);
areValid = areValid & (length(Rfi) >= itNum);
areValid = areValid & (length(Rsi) >= itNum);
areValid = areValid & (length(Rsai) >= itNum);
end % valuesAreValid function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hasComplex] = checkForComplexGoalResType()
hasComplex = false;
for count = 1:Nr
    % TODO_DWW: Use regex split to just grab {"_complex"} generalising to any s-param
    if strcmp(OPTopts.goalResType{count},'S11_complex')
        hasComplex = true;
        break
    end
end % for Nr
end % valuesAreValid function

end % plotModels function