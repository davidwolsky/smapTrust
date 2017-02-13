
function plotModels(plotFlag, itNum, Rci, Rfi, Rsi, Rsai, OPTopts)
    if plotFlag && valuesAreValid()
        figure(itNum)
		Nr = length(OPTopts.Rtype); % Number of responses requested
        for rr = 1:Nr
            subplot(Nr,1,rr)
            if isfield(Rci{itNum}{rr},'f')
                plot(Rci{itNum}{rr}.f,Rfi{itNum}{rr}.r,'k'), grid on, hold on
                plot(Rci{itNum}{rr}.f,Rci{itNum}{rr}.r,'r'), grid on, hold on
                plot(Rsi{itNum}{rr}.f,Rsi{itNum}{rr}.r,'b')
                plot(Rsai{itNum}{rr}.f,Rsai{itNum}{rr}.r,'g--')
                xlabel('Frequency')
                % Plot the specs...
            else
                plot(Rfi{itNum}{rr}.r,'k'), grid on, hold on
                plot(Rci{itNum}{rr}.r,'r'), grid on, hold on
                plot(Rsi{itNum}{rr}.r,'b')
                plot(Rsai{itNum}{rr}.r,'g')
                xlabel('Index')
                % Plot the specs...
            end
            ylabel(OPTopts.Rtype{rr})
            title(['Iteration ',num2str(itNum)])
            legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate')
        end

    end
    
    function [areValid] = valuesAreValid()
        areValid = true;
        areValid = areValid & (length(Rci) >= itNum);
        areValid = areValid & (length(Rfi) >= itNum);
        areValid = areValid & (length(Rsi) >= itNum);
        areValid = areValid & (length(Rsai) >= itNum);
    end % 

end % plotModels











