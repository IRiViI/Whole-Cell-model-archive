% Template for making new function for the whole cell model archive class.
%   ARCHIVE = TEMPLATE_ARCHIVE(ARCHIVE) run simulation with default 
%   settings.
%
%   'set', sets                 - Specify a specific set of sets to analyze
%
%   'simulation', simulations   - Specify specific simulations to analyze
%
%   'progress', state           - Turn on or off progress report (true or
%                                 false)
%
%   'file', file                - Change the file to open in one of the
%                                 examples
%
%   ARCHIVE = TEMPLATE_ARCHIVE(ARCHIVE, 'progress', FALSE) no progress
%   report.
%
%   ARCHIVE = TEMPLATE_ARCHIVE(ARCHIVE,'set',SETS,'simulation',SIMULATIONS)
%   Only process a specific set and/or simulation
%
%   Notes:
%
%   There are two areas highlighted in this document which are designated
%   for adding your own code. These areas contain examples which can be
%   removed without consequences.
%
%   The minimal number of files required in the simulation folders:
%   options.mat, summary.mat, state-0.mat and state-1.mat.
%
%   Directory recommanded:
%   WholeCell-master
%
%   This function is functional after executing the following procedure:
%   >>  archive = Archive();        % Create new archive object.
%   >>  archive.add_sets(path,'short',true);
%                                   % path: path to a directory that
%                                   contains folders with simulation files.
%                                   example: '/home/path/to/dir'
%                                   "'short',true" optional but recommended
%   >>  archive.initiate_WCM;       % Recommended, but not necessary.
%   >>  archive.compress_state_files 
%                                   % Compress the state files.
%   >>  archive.general_processing  % Extract some information of the state
%                                   files.
%   Finally:
%   >>  archive.template
%   
%   This function can be implemented in the class "Archive.m" as follows:
%   methods
%       function this = template(this,varargin)
%           this = template_archive(this,varargin{:});
%       end
%   end
%
%   After completion of a new method, place the new script in the
%   subfolder "methods"
%
%   Examples:
%
%   template_archive(archive, 'set', 1, 'progress',false)
%       Only analyze the first set of the archive structure without giving
%       a progress report
%
%   [archive, example] = template_archive(archive, 'file', 'state-2.mat');
%       Analyze a different file in on of the examples. Get the example
%       output as well as the archive output.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function varargout = template_archive(archive, varargin)
%% Settings

% Default settings
options.set = 0;        % I usually use a value of 0 to indicate that all sets should be analyzed
options.simulation = 0; % Analyze all the simulations of the sets.
options.progress = true;    % Show progress

% Default setting for example:
options.file = 'point_compressed_state.mat';

% Change settings
tmpOptions = struct(varargin{1:end});                       % Make a structure of the varargin inputs
field = fields(tmpOptions);                                 % Get all the fiels of teh varargin structure
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});   % Set all the values of the input structure to the existing options structure
end

% Note:
% archive is the first input on default. It's possible to add other default 
% inputs like: 
% template_archive(archive, input2, input3)
% In this case adjust the function accordingly:
% varargout = template_archive(archive, input2, input3, varargin)

%% Process settings

% Determine which sets of the archive should analyzed
if options.set == 0
    tSet = archive.sets;    % Get the total number of sets
    lSet = 1:tSet;          % List of sets that will be analyzed includes all sets
else
    lSet = options.set;     % List of sets that will be analyzed includes all the sets according to the input
    tSet = length(lSet);    % Get the total number of sets that will be analyzed in this function. Note: THIS DOES NOT HAVE TO BE THE SAME NUMBER OF SIMULATIONS THAT EXISTS IN THE ARCHIVE
end

% Determine the total number of simulations that will be analyzed in this
% function (this is used for showing the progress of the function)
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

% For every set that whould be analyzed
for iSet = lSet
    
    % Determine which simulations of the set should analyzed
    if options.simulation == 0
        tSimulation = archive.simulations(iSet);    % Get the total number of simulations of the current set
        lSimulation = 1:tSimulation;                % List of simulations that should be analyzed includes all the simulations of the current set
    else
        lSimulation = options.simulation;                  % List of simulations that will be analyzed includes all the simulations according to the inputs
        tSimulation = length(lSimulation);          % Number of simulations to analyze
    end
    
    
    % ---------------------------------------------------------------%
    % -------- YOUR CODE FOR ANALYZING INDIVIDUAL SETS HERE -------- %
    % \/----------\/----------\/----------\/----------\/----------\/ %
    
    
    
    
    
    % Examples: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ begin
    
    % Get type of knockout
    knockoutGene = archive.set(iSet).knockout;
    
    % Get trend values of the set
    massTrendSet = [archive.set(iSet).trends.average_time];
    massTrendSet = [archive.set(iSet).trends.average_mass];
    
    % Set new values, example 1:
    %     archive.set(iSet).newValue = 'newValue';
    
    % Set new values, example 2:
    newValue = struct;
    newValue.field1 = 'newValue';
    newValue.field2 = [1,2,3];
    newValue(2).field1 = 'newValue';
    newValue(2).field2 = [4,5,6];
    %     archive.set(iSet).newValue = newValue;
    
    % Examples: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end
    
    
    
    
    
    % /\----------/\----------/\----------/\----------/\----------/\ %
    % -------- YOUR CODE FOR ANALYZING INDIVIDUAL SETS HERE -------- %
    % -------------------------------------------------------------- %
    
    
    % For every simulation in the set
    for iSimulation = lSimulation
        
        
        % ---------------------------------------------------------------%
        % ---- YOUR CODE FOR ANALYZING INDIVIDUAL SIMULATIONS HERE ----- %
        % \/----------\/----------\/----------\/----------\/----------\/ %
        
        
        
        
        
        
        % Examples: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ begin
        
        % Get folder of current simulation:
        folder = archive.set(iSet).simulation(iSimulation).folder;
        
        % Load the fields "Time" and "ProteinComplex" of a specific state
        % file
        number = 1;
        file = ['state-' num2str(number) '.mat'];
        field = {'Time','ProteinComplex'};
        state = archive.load_file(iSet,iSimulation,...
            'fileTypeTag',file,...  % Optional (default: load point_compressed_state.mat)
            'fields',field);         % Optional (default: load all the fields of the file)
        
        % Load compressed state file according to the options.
        file = options.file;
        state = archive.load_file(iSet,iSimulation,...
            'fileTypeTag',file);   % Optional (default: load point_compressed_state.mat), this option is redundant as well.
        
        % Load options file (be careful not to overwrite the variable
        % "options")
        file = 'options.mat';
        optionFile = archive.load_file(iSet,iSimulation,...
            'fileTypeTag',file);   % Optional (default: load point_compressed_state.mat)
        
        % Access a value of the archive structure
        massTrend = [archive.set(iSet).simulation(iSimulation).trends.mass];
        timeTrend = [archive.set(iSet).simulation(iSimulation).trends.time];
        seedSim = archive.set(iSet).simulation(iSimulation).seed;
        
        % Set new values, example 1:
        %     archive.set(iSet).simulation(iSimulation).newValue = 'newValue';
        
        % Set new values, example 2:
        newValue = struct;
        newValue.field1 = 'newValue';
        newValue.field2 = [1,2,3];
        newValue(2).field1 = 'newValue';
        newValue(2).field2 = [4,5,6];
        %     archive.set(iSet).simulation(iSimulation).newValue = newValue;
        
        % Examples: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end
        
        
        
        % /\----------/\----------/\----------/\----------/\----------/\ %
        % ---- YOUR CODE FOR ANALYZING INDIVIDUAL SIMULATIONS HERE ----- %
        % -------------------------------------------------------------- %
        
        
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

% Set outputs
varargout{1} = archive;
varargout{2} = 'example';

end

