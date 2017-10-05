function unitString = extractUnitString(goalResType)

% Function to extract a unit type from the goalResType from OPTopts. The 
% unit string is separated uniquely by an underscore '_'. 
% The unit string is returned.

% TODO_DWW: See if find(...) works better here else expand.
splits = regexp(goalResType, '\_', 'split');
assert(length(splits) == 2);
unitString = splits{2};

end % extractUnitString function