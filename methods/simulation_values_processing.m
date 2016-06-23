% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/27/2016

function varargout = simulation_values_processing(archive, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.fileTypeTag = 'state-';
options.target = {'Metabolite.counts',...
    'ProteinComplex.counts',...
    'ProteinMonomer.counts',...
    'RNA.counts',...
    'MetabolicReaction.fluxs'};
options.timeStep = 1;

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

% Number of targets
tTarget = size(options.target,2);

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
        
        % Process the targets
        for iTarget = 1:tTarget
           
            % Current target
            target = options.target{iTarget};
            
            if strcmp(target,'MetabolicReaction')
                compartment = 1;
            else
                compartment = 1:6;
            end
            
            offset = 50;
            eRange = 10;
            
            timeStart = 10000;
            fileInfoStart = archive.find_time(timeStart);
            timeStop = 10999;
            fileInfoStop = archive.find_time(timeStop);
            
            elements = 1+offset:eRange+offset+1;
            
            matrix.time = fileInfoStart.Time.min:options.timeStep:fileInfoStop.Time.max;
            
            % Get the values of a specific quantity over the whole
            % simulation.
            matrix.values = archive.get_values(iSet,iSimulation,target,elements,compartment,...
                'stateStart',fileInfoStart.number,...
                'stateStop',fileInfoStop.number);
            
            
            element = 54;
            iElement = element - offset;
            location = 1;
            
            time = 100;
            range = 3;
            
            values = [time-range:time+range];
            
            iMatrix = matrix(iElement,location,values);
            
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
varargout{1} = archive;

end