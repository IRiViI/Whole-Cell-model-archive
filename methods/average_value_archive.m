% [strainAverage, simulationValues] = average_value_archive(archive,...
%   strain,timePoint,field,subfield,location);
%
%   'fileTypeTag'   - Give a unique part of the file which should be used
%                     to obtain the value. default: 'mean_compressed'.
%   'progress'      - Show progress of getting value. Default: 'Off' 
%

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/05/2016

function [strainAverage,matrix] = average_value_archive(varargin)

% Load optional inputs
if nargin > 6
   options = struct(varargin{7:end}); 
end

% Process input arguments
archive = varargin{1};
options.strain = varargin{2};
options.timePoint = varargin{3};
options.field = varargin{4};
options.subfield = varargin{5};
options.location = varargin{6};

% Add default settings if no values are given
if ~isfield(options,'fileTypeTag')
    options.fileTypeTag = 'mean_compressed_state.mat';
end
if ~isfield(options,'progress')
    options.progress = 'Off';
end
if ~isfield(options,'direct')
    options.direct = 'Off';
end

% Number of simulations in strain
tSimulation = length(archive.set(options.strain).simulation);

% Progress display offset
if strcmp(options.progress,'On')
    fprintf('Progress:\n     ');
end

% Number of time points
tTime = length(options.timePoint);

% Number of locations
t_location = length(options.location);

% For every simulation
for iSimulation = 1:tSimulation
    
    % Display progress
    if strcmp(options.progress,'On')
        display_progress(iSimulation,tSimulation)
    end
    
    % Load state file
    if strcmp(options.direct,'On')
        state = load([archive.set(options.strain).simulation(iSimulation).folder '/' options.fileTypeTag],options.field,'Time');
    else
        state = load_state(archive,options.strain,iSimulation,...
            'fileTypeTag',options.fileTypeTag,...
            'fields',{options.field,'Time'});
    end
    
    for iTime = 1:tTime
        % Get time points
        [~,timePoint] = min(abs(ceil(state.Time.values - options.timePoint(iTime))));
        for i_location = 1:t_location
            location = options.location(i_location);
        % Load values to one big matrix
        matrix(:,i_location,iTime,iSimulation) = squeeze(state.(options.field).(options.subfield)(:,location,timePoint));
    
        end
    end
end

% Get average value
for iTime = 1:tTime
    strainAverage(:,iTime,:) = mean(matrix(:,:,iTime,:),4);
end
end
