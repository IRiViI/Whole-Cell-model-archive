% Template for making new function for the whole cell model archive class.
%   ARCHIVE = TEMPLATE_ARCHIVE(ARCHIVE) run simulation with default 
%   settings.
%   Stripped version of "template_archive.m".

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function varargout = empty_template_archive(archive, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress

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
        lSimulation = options.set;                  % List of simulations to analyze
        tSimulation = length(lSimulation);          % Number of simulations to analyze
    end
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Code here (analyzes of different sets)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Analyze simulations
    for iSimulation = lSimulation
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Code here (analyzes of different simulations
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
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
varargout{1} = archive;

end

