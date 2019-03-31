function populateSpaceAndPlot(prepopulatedSpaceFile, plotModelOpts, plotIterationOpts)

% Takes a matlab data file generated by the logging system (createLog), loads it and 
% uses the plotting functions to regenerate the plot. Basic images are generated for the
% figures.
 
% Inputs:
%   prepopulatedSpaceFile: File to load the varaible data space from.
%   plotGoalFlag: If false then the goal values are not added to the plot.
%   plotModelOpts: Options for the model plotting.
%   plotIterationOpts: Options for the iteration plotting.

% TODO:
%   - Add figures for the alignment phase.

plotIter = 1;
plotGoalFlag = 0;

if isfield(plotModelOpts, 'plotGoalFlag'), plotGoalFlag = plotModelOpts.plotGoalFlag; end

space = getPrepopulatedSpaceFrom(prepopulatedSpaceFile);

% keyboard
% TODO_DWW: Check how it can be that this sometimes has one less ii item than its meant to have.

for count = 1:1:space.ii
    plotModels(plotIter, plotGoalFlag, count, space.Rci, space.Rfi, space.Rsi, ...
            space.Rsai, space.OPTopts, plotModelOpts);
    print(['plotModel_', mat2str(count) ], '-depsc')
end

plotNormalised = 1;
plotIterations(plotIter, space.xi, space.Ti.Delta, space.OPTopts, space.SMopts, ...
        space.Si, plotNormalised, 'Normalised', plotIterationOpts);
print(['plotIteration_normalised', mat2str(count) ], '-depsc')
plotIterations(plotIter, space.xi, space.Ti.Delta, space.OPTopts, space.SMopts, space.Si, ...
         ~plotNormalised, 'De-normalised/globalised/universalised', plotIterationOpts);
print(['plotIteration_deNormalised', mat2str(count) ], '-depsc')

end % populateSpaceAndPlot