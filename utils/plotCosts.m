function plotCosts(Ti, OPTopts, costS, costF)

% Plot the costs (typically after SM run) for the fine and 
% surrogate models as functions of the number of iterations.

% Inputs:
%   Ti:         Trust region structure.
%   OPTopts:    Optimization options, used for determining what the type of
%               response is and goal values, as well as get all the costs used in the
%               TR
%   costS:      Surrogate model costs per iteration.
%   costF:      Fine model costs per iteration.

markerstr = 'xso+*d^v><ph.xso+*d^v><ph.';
colourstr = 'kbrgmcykbrgmckbrgmcykbrgmc';

figure()
Ni = length(Ti.xi_all);    % Number of iterations
% costSi1 = [];
% for nn = 1:Ni
%     for cc = 1:Ni
%         costSi1(cc) = costSurr(Ti.xin_all{cc},Ti.Si_all{nn}{:},OPTopts);
%     end
%     subplot(2,2,1)
%     plot(costSi1, strcat(markerstr(nn),colourstr(nn)),'LineWidth',2,'MarkerSize',10), grid on, hold on
%     ylabel('costSi')
%     xlabel('Iterations all')
% end
% legend('show') 
% title('costs')
subplot(2,2,[1 2])
if strcmp(OPTopts.goalResType{1},'Gen')
    plot(cell2mat(Ti.costF_all), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'complex')
    plot(real(cell2mat(Ti.costF_all)), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plot(imag(cell2mat(Ti.costF_all)), strcat(markerstr(2),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    legend('real','imag');
    plotRealGoals()
    plotImagGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'dB')
    plot(cell2mat(Ti.costF_all), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
else
    error(['goalResType (', OPTopts.goalResType{1}, ')not found']);
end
ylabel('Ti.costF\_all')
xlabel('Iterations all')
title('Ti.costF\_all')

% Meaningless -> just for debugging
subplot(2,2,3)
if strcmp(OPTopts.goalResType{1},'Gen')
    plot(cell2mat(costS), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'complex')
    plot(real(cell2mat(costS)), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plot(imag(cell2mat(costS)), strcat(markerstr(2),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    legend('real','imag')
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'dB')
    plot(cell2mat(costS), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
else
    error(['goalResType (', OPTopts.goalResType{1}, ')not found']);
end
ylabel('costS')
xlabel('Iterations success points')
title('costS - fairly meaningless')

subplot(2,2,4)
if strcmp(OPTopts.goalResType{1},'Gen')
    plot(cell2mat(costF), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'complex')
    plot(real(cell2mat(costF)), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plot(imag(cell2mat(costF)), strcat(markerstr(2),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    legend('real','imag')
    plotRealGoals()
    plotImagGoals()
elseif strcmp(extractUnitString(OPTopts.goalResType{1}),'dB')
    plot(cell2mat(costF), strcat(markerstr(1),colourstr(1)),'LineWidth',2,'MarkerSize',10), grid on, hold on
    plotGoals()
else
    error(['goalResType (', OPTopts.goalResType{1}, ')not found']);
end
ylabel('costF')
xlabel('Iterations success points')
title('costF')


% ----- goals -----
function plotGoals()
if ( isfield(OPTopts, 'goalVal') )%&& isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        plot([1, Ni], (OPTopts.goalVal{gg})*ones(1,2), 'm', 'LineWidth',2)
    end
end % if validation
end % plotGoals function

function plotRealGoals()
if ( isfield(OPTopts, 'goalVal') )%&& isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        plot([1, Ni], real(OPTopts.goalVal{gg})*ones(1,2), 'm', 'LineWidth',2)
    end
end % if validation
end % plotRealGoals function

function plotImagGoals()
if ( isfield(OPTopts, 'goalVal') )%&& isfield(OPTopts, 'goalStart') && isfield(OPTopts, 'goalStop') )
    Ng = length(OPTopts.goalVal);
    for gg = 1:Ng
        if ( ~isreal(OPTopts.goalVal{gg}) )
            plot([1, Ni], imag(OPTopts.goalVal{gg})*ones(1,2), 'c', 'LineWidth',2)
        end % imag part
    end
end % if validation
end % plotImagGoals function
% -----  -----

end % plotCosts