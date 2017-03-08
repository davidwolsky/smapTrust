
function plotModels(plotFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts)
if plotFlag && valuesAreValid()
    figure(itNum)
	Nr = length(OPTopts.Rtype); % Number of responses requested
    for rr = 1:Nr
        subplot(Nr,1,rr)
        if isfield(Rci{itNum}{rr},'f')
            plot(Rci{itNum}{rr}.f,Rfi{itNum}{rr}.r,'k','LineWidth',1.5), grid on, hold on
            plot(Rci{itNum}{rr}.f,Rci{itNum}{rr}.r,'r','LineWidth',1.5), grid on, hold on
            plot(Rsi{itNum}{rr}.f,Rsi{itNum}{rr}.r,'b','LineWidth',1.5), grid on, hold on
            plot(Rsai{itNum}{rr}.f,Rsai{itNum}{rr}.r,'g--','LineWidth',2), grid on, hold on
            xlabel('Frequency')
        else
            plot(Rfi{itNum}{rr}.r,'k','LineWidth',1.5), grid on, hold on
            plot(Rci{itNum}{rr}.r,'r','LineWidth',1.5), grid on, hold on
            plot(Rsi{itNum}{rr}.r,'b','LineWidth',1.5), grid on, hold on
            plot(Rsai{itNum}{rr}.r,'g','LineWidth',1.5), grid on, hold on
            xlabel('Index')
        end
        plotGoals()
        ylabel(OPTopts.Rtype{rr})
        title(['Iteration ',num2str(itNum)])
        legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate')
    end
end % if validation


function plotGoals()
if ( isfield(OPTopts, 'goalVal') && isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        plot([OPTopts.goalStart{gg}, OPTopts.goalStop{gg}], OPTopts.goalVal{gg}*ones(1,2), 'm', 'LineWidth',3)
    end
end % if validation
end % plotGoals

function [areValid] = valuesAreValid()
areValid = true;
areValid = areValid & (length(Rci) >= itNum);
areValid = areValid & (length(Rfi) >= itNum);
areValid = areValid & (length(Rsi) >= itNum);
areValid = areValid & (length(Rsai) >= itNum);
end % valuesAreValid

end % plotModels