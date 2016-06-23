function varargout = snap_shot_average_archive(archive, snap_shot, sets)
% snap_shot_average_archive(archive, snap_shot,sets)
%   [meanMatrix, stdMatrix] = snap_shot_average_archive(archive,
%   snap_shot, sets)

% Initiate structure
structure = struct([]);

% Process all the sets of a set
for iSet = sets
    
    % Determine the total number of simulations in set
    tSimulation = archive.simulations(iSet);
    
    % Add all values of the simulations to structure
    for iSimulation = 1:tSimulation
        structure(end+1).values = archive.set(iSet).simulation(iSimulation).snapShot(snap_shot).values;
    end
    
end

% Get number of dimension
dimensions = length(size(structure(1).values));

% Create a single matrix
total_structure = cat(dimensions+1,structure.values);

% Determine average and std
varargout{1} = mean(total_structure,dimensions+1);
varargout{2} = std(total_structure,[],dimensions+1);

end