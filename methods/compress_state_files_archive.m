% Compress the output files of the whole cell model. This script combines
%   the data of different 'state' files in two different ways. 1) last data
%   points of all the 'state' files of a folder are saved. 2) The mean
%   value is calculated of all the values in a single state file.
%   In this version are multiple simulation folders processed. Select
%   the directory that contains folders which contain the state files.
%
%   compress_state_files_archive(archive)
%
%   'set'           - Select specific set. Set set to "0" to process all
%                     sets
%   'simulation'    - Select specific simualtion. Set simulation to "0" to
%                     process all simulations
%   'redo'          - redo simultions (true or false). default: false
%   'core'          - UNDER CONSTRUCTION
%   'minimal'       - Compress simulations without using the special
%                     libraries of the whole cell model. Note: This will
%                     result is a less complete point_compressed_state file
%   'warning'       - Turn off warnings ('warning',false)
%   'externalMultiCore'
%                   - Only compress one simulation specific simulation of 
%                     the compression library. This is usefull when
%                     compressing files on a cluster.
%
%   Examples:
%       SBATCH ARRAY                                                      
%       #SBATCH --array=1-100    
%       compress_state_files_archive('externalMultiCore',n)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 10/05/2016

function compress_state_files_archive(varargin)

%% Process input arguments

% Mandatory input
archive = varargin{1};

% Default settings
options.set = 0;
options.simulation = 0;
options.core = 1;
options.redo = false;
options.warning = true;
options.minimal = false;
options.externalMultiCore = 0;

% Adjust options
inputOptions = struct(varargin{2:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

%% Process settings

% Prepare parallel processing
if options.core > 1
    % Get informatoin current parpool
    p = gcp('nocreate');
    % Check if parpool should be created/adjusted
    if isempty(p) % If parpool does not exist
        p = parpool(options.core);
    elseif p.NumWorkers ~= options.core % If parpool number of workes is not correct
        delete(p);
        p = parpool(options.core);
    end
end

%% Compressing

% Sets to process
tSet = archive.sets;
if options.set == 0
    lSet = 1:tSet;
else
    lSet = options.set;
end

% Total number of simulations processed
if options.simulation == 0
    ttSimulation = archive.simulations(lSet);
else
    ttSimulation = length(lSet) * length(options.simulation);
end

% Progress
iiSimulation = 0;
fprintf('Progress:\n     ')
display_progress(iiSimulation,ttSimulation);

% Initate simulation library (This makes it easier to to use muliple cores)
libSimulation = [];

% Generate simulation library
for iSet = lSet
    
    % Simulations to process
    tSimulation = archive.simulations(iSet);
    if options.simulation == 0
        lSimulation = 1:tSimulation;
    else
        lSimulation = options.simulation;
    end
    
    % Put all the targeted simulations of the targeted set in the library
    for iSimulation = lSimulation
        
        % Generate library for all the simulation targets
        libSimulation(end+1,:) = [iSet,iSimulation];
        
    end
    
end

% All the simulations targeted
tSim = length(libSimulation);

if options.core > 1 % If multicore
    
    % Turn off warning messages
    if ~options.warning
        pctRunOnAll warning off
    end
    
    % Process all the simulations in the library
    parfor iSim = 1:tSim
        
        % Get the simulation and set number
        iSet = libSimulation(iSim,1);
        iSimulation = libSimulation(iSim,2);
        
        % Folder
        folder = archive.set(iSet).simulation(iSimulation).folder;
        
        % Determine state files
        [endState,meanState] = create_structures_from_files_all_fields(folder);
        
        % Save
        parsave([folder '/' 'mean_compressed_state.mat'],'-struct','meanState');
        parsave([folder '/' 'point_compressed_state.mat'],'-struct','endState');
        
        fprintf('Done: set %d simulation %d (folder %s)\n',iSet,iSim,folder);
    end
    
else % If not multicore
    
    % Turn off the warnings
    if ~options.warning
        warning off;
    end
    
    % Range of library to compress
    lSim = 1:tSim;
    
    % Select only one simultion of the library when dealing with a mulit
    % core processor which is not external of the matlab program
    if options.externalMultiCore ~= 0
        lSim = options.externalMultiCore;
    end
    
    % Process all the simulations in the library
    for iSim = lSim
        
        % Get the simulation and set number
        iSet = libSimulation(iSim,1);
        iSimulation = libSimulation(iSim,2);
        
        % Folder
        folder = archive.set(iSet).simulation(iSimulation).folder;
        
        % Determine state files
        [endState,meanState] = create_structures_from_files_all_fields(folder);
        
        % Save
        save([folder '/' 'mean_compressed_state.mat'],'-struct','meanState');
        save([folder '/' 'point_compressed_state.mat'],'-struct','endState');
        
        % Progress
        iiSimulation = iiSimulation + 1;
        display_progress(iiSimulation,ttSimulation);
    end
    
end

fprintf('Done\n');

end

function [endState,meanState] = create_structures_from_files_all_fields(folder,options)
% Get the all the compressed states of the folder

iTime = 1; % Frame counter

% Check every file of folder
while true
    
    % Update state file name
    nState = ['state-' num2str(iTime - 1) '.mat'];
    
    % load next state file
    try
        state = load([folder '/' nState]);
    catch
        % State file does not exist
        break
    end
    
    % Get field information
    fieldNames = fieldnames(state);
    tField = length(fieldNames);
    
    % Get all the data every field
    for iField = 1:tField
        
        % Get sub-field information
        subFieldNames = fieldnames(state.(fieldNames{iField}));
        tSubField = length(subFieldNames);
        
        % For every subfield
        for iSubField = 1:tSubField
            try
                % Method 1: Process normal structures
                % Get the last value of every data type
                try
                    endState.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,iTime) = state.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,end);
                catch
                    % This method is sometimes require when
                    % there aren't 3 dimensions
                    endState.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,iTime) = state.(fieldNames{iField}).(subFieldNames{iSubField});
                end
                meanState.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,iTime) = mean(state.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,:),3);
            catch
                if ~options.minimal % If all the information should be included
                    try
                        % Method 2: Process special structures that
                        % contain separate fields
                        try
                            endState.(fieldNames{iField}).(subFieldNames{iSubField})(iTime) = struct(state.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,end));
                        catch
                            endState.(fieldNames{iField}).(subFieldNames{iSubField})(iTime) = struct(state.(fieldNames{iField}).(subFieldNames{iSubField}));
                        end
                    catch
                        % Method 3: Process special structures that
                        % do not contain serpate fields
                        try
                            endState.(fieldNames{iField}).(subFieldNames{iSubField})(iTime).data = state.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,end);
                        catch
                            try
                                endState.(['extra_' fieldNames{iField}]).(subFieldNames{iSubField})(:,:,iTime) = state.(fieldNames{iField}).(subFieldNames{iSubField})(:,:,end);
                            catch
                                warningMessage = sprintf('Unable to process field %s subfield %s of file %s folder %s',...
                                    (fieldNames{iField}), (subFieldNames{iSubField}), nState, folder);
                                warning(warningMessage);
                            end
                        end
                    end
                end
            end
            
        end
        
    end
    
    iTime = iTime + 1;
    
end

% varargout{1} = endState;
% varargout{2} = meanState;

end

function parsave(varargin)

% save to change the '-struct' components of saving or something...

save(varargin{:});
end