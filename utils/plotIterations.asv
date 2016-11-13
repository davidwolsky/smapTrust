function plotIterations(plotFlag, xi, Delta, OPTopts, plotNormalised, caption)
% Plotter for iteration by iteration details. 
% For a two dimensional problem the iterations will be plotted against each other. For more than two 
% parameters the plots take place iteration by iteration. The trust region radius (assuming square shape)
% is plotted for each iteration.

% Inputs:
%   plotFlag: 
%   xi: The parameters.
%   Delta: The trust region radius at each iteration. 
%   OPTopts:  - ximin and ximax are plotted. 
%             - Pass a 'false' if these are not required.
%   plotNormalised: use OPTopts ximin & ximax to normalise the parameters to plot. 
%   caption: The caption of the graph.

if plotFlag
    markerstr = 'xso+*d.^v><ph';
    colourstr = 'kbrgmcy';

    % Manipulate data format and ensure that everything has the same dimensions
    Ni = length(xi);    % Number of iterations
    Nx = length(xi{1}); % Number of parameters
    
    Ndi = length(Delta);    %Number of TR radius entries
    Ndx = length(Delta{1}); %Number of TR radius parameters
    if ( Ni > Ndi )
        % Handle the case of the no TR is available - returned from tolerance or something 
        % TODO_DWW: Check if this is fixed
        Delta{end+1} = zeros(Ndx,1);
    end
    optHasBounds = ( isfield(OPTopts, 'ximax') && isfield(OPTopts, 'ximax') );
    if plotNormalised && optHasBounds
        for ii = 1:Ni
            xi{ii} = (xi{ii} - OPTopts.ximin)./(OPTopts.ximax - OPTopts.ximin);
            Delta{ii} = (Delta{ii})./(OPTopts.ximax - OPTopts.ximin);
        end
    end
    transDelta = transpose(cell2mat(Delta));
    transXi = transpose(cell2mat(xi));

    figure()
    % For two dimensional plot the two parameters against each other 
    if ( Nx == 2 )
        axis equal
        for ii = 1:Ni
            plot(xi{ii}(1),xi{ii}(2),strcat(markerstr(1),colourstr(ii)),'LineWidth',2,'MarkerSize',5*ii), grid on, hold on
            % Plot the TR radius by setting up the rectangle
            % if isfield(OPTopts, 'TREnabled') && OPTopts.TREnabled
                if (Ndx > 1)
                    transDDelta = [transDelta];
                else
                    transDDelta = [transDelta,transDelta];
                end
                transMerge = [transXi-transDDelta,transDDelta*2];
                % Only plot radius if it is valid
                if (transMerge(ii,3) > 0)
                    rectangle('Position',transMerge(ii,1:4), 'EdgeColor',colourstr(ii),'LineWidth',2)
                end

                % Plot optimiser min and max
                if optHasBounds
                    if plotNormalised
                        rectangle('Position',[0,0,1,1])
                    else
                        for ii = 1:Nx
                            rectangle('Position',[transpose(OPTopts.ximin),transpose(abs((OPTopts.ximin))+OPTopts.ximax)])
                        end
                    end
                end

            % end
        end
        axis equal
        xlabel('x1 value')
        ylabel('x2 value')
        title({strcat('Values at each iteration plotted against each other with trust region radius - ', caption)})
    else
        for ii = 1:Nx
            plot(transXi(:,ii), strcat(markerstr(ii),colourstr(ii)),'LineWidth',2,'MarkerSize',10), grid on, hold on
            % if isfield(OPTopts, 'TREnabled') && OPTopts.TREnabled
                xOffset = transpose(1:Ni)+(ii*0.07);
                if (Ndx > 1)
                    errorbar(xOffset, transXi(:,ii), transDelta(:,ii), strcat('.',colourstr(ii)),'LineWidth',1.5,'MarkerSize',10 ), grid on, hold on
                else
                    errorbar(xOffset, transXi(:,ii), transDelta(:), strcat('.',colourstr(ii)),'LineWidth',1.5,'MarkerSize',10 ), grid on, hold on
                end
            % end
        end
        xlabel('Iteration')
        ylabel('Value')
        title({strcat('Values per iteration with trust region radius - ', caption)})
        
        % Plot optimiser min and max
        if optHasBounds
            if plotNormalised
                    plot(ones(1,Ni))
                    plot(zeros(1,Ni))
            else
                for ii = 1:Nx
                    plot(OPTopts.ximax(ii)*ones(1,Ni))
                    plot(OPTopts.ximin(ii)*ones(1,Ni))
                end
            end
        end
    end
    legend()
end
end % plotIterations
