% Extract data the elements of different simulations. 

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function varargout = sim_set_matrix_archive(archive, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.element = 0;       % Mandatory
options.field = '';         % Mandatory input
options.subfield = '';      % Mandatory input
options.minTime = 1;
options.maxTime = inf;
options.compartment = inf;
options.fileTypeTag = 'mean_compressed_state.mat';

% Change settings
tmpOptions = struct(varargin{1:end});                      
field = fields(tmpOptions);                                 
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});   
end

%% Process settings

% Determine which sets of the archive should be analyzed
if options.set == 0
    tSet = archive.sets;    % Number of sets to analyze
    lSet = 1:tSet;          % List of sets to analyze
else
    lSet = options.set;     % List of sets to analyze
    tSet = length(lSet);    % Number of sets to analyze
end

% Determine number of simulations to analyze (required for progess report)
if options.simulation == 0
    ttSimulation = archive.simulations(lSet);
else
    ttSimulation = tSet * length(options.simulation);
end

% Check if all mandatory inputs are given
if ~isfield(options,'field') || ...
        ~isfield(options,'subfield') ||...
        ~isfield(options,'element')
   error('Mandatory inputs includes: field, subfield and elements') 
end

% Number of compartments when not defined
if options.compartment == inf
    folder = archive.folder(lSet(1),1);
    state = load([folder '/' options.fileTypeTag],options.field);
    values = state.(options.field).(options.subfield);
    % Include all compartments in calculations
    options.compartment = 1:size(values,2);
end

%% Run code

% Initiate simulation counter
iiSimulation = 0; % Simulation counter

% Display progress
if options.progress
    fprintf('Progress:\n     ');
    display_progress(iiSimulation, ttSimulation)
end

% Analyze sets
for iSet = lSet
    
    % Make list of sets to analyze
    if options.simulation == 0
        tSimulation = archive.simulations(iSet);    % Number of simulations to analyze
        lSimulation = 1:tSimulation;                % List of simulations to analyze
    else
        lSimulation = options.simulation;                  % List of simulations to analyze
        tSimulation = length(lSimulation);          % Number of simulations to analyze
    end
    
    % Analyze simulations
    for iSimulation = lSimulation
        
        % Get folder
        folder = archive.folder(iSet,iSimulation);
        
        % Load State file
        state = load([folder '/' options.fileTypeTag],options.field);
        
        % Extract values
        if options.maxTime == inf
            values = state.(options.field).(options.subfield)(options.element,options.compartment,options.minTime:end);
        else
            values = state.(options.field).(options.subfield)(options.element,options.compartment,options.minTime:options.maxTime);
        end
        
        % Add values to output structure
        if exist('totalValues','var')
            totalValues = cat(4,totalValues,values);
        else
            totalValues = values;
        end
        
        % Update simulation counter
        iiSimulation = iiSimulation + 1;
        
        % Display progress
        if options.progress
            display_progress(iiSimulation, ttSimulation)
        end
        
    end
    
end

% Display progress
if options.progress
    fprintf('Done\n');
end

%% Set outputs

% Set output
varargout{1} = totalValues;
varargout{2} = options;

end

