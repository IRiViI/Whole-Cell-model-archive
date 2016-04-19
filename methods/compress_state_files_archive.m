% Compress the output files of the whole cell model. This script combines 
%   the data of different 'state' files in two different ways. 1) last data 
%   points of all the 'state' files of a folder are saved. 2) The mean 
%   value is calculated of all the values in a single state file. 
%   In this version are multiple simulation folders processed. Select
%   the directory that contains folders which contain the state files.
%
%   Example:
%       compress_state_files_archive('dataDir','path/to/simulations/folder')
%
%   Other example:
%       compress_state_files_archive('dataDir','path/to/simulations/folder',...
%           'outputDir','path/to/output/folder',...
%           'tag','state',...
%           'outputName','begin_of_name_of_output_file',...
%           'data',A,...
%           'subFolder','Off')
%               A.name = {'Rna','Time','Metabolite','ProteinComplex',...
%                           'ProteinMonomer','Mass'};

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/04/2016

function compress_state_files_archive(varargin)

%% Process input arguments

% Get optional inputs
if nargin > 1
    options = struct(varargin{2:end});
else
    options = struct();
end
% Mandatory input
options.dataDir = varargin{1};

% Set output dir if not given
if ~isfield(options,'outputDir')
    options.outputDir = options.dataDir; 
end

% The name of the files that contain the values of the simulations
if ~isfield(options, 'tag')
    % Default
    options.tag = 'state-';
end

% The fields of the structures that will be transfered
if ~isfield(options, 'data')
    % Default
    %     options.data.name = {'Rna','Time','Metabolite','ProteinComplex','ProteinMonomer','Mass'};
    options.data.name = {'All'};
end

% First part of the name of the output file
if ~isfield(options, 'outputName')
    % Default
    options.outputName = 'compressed_state';
end

% Set whetehr if the subfolders should be conserved or not
if ~isfield(options, 'subFolder')
    % Default
    options.subFolder = 'On';
end

% Compress the files of a single folder
if ~isfield(options, 'singleFolder')
    options.singleFolder = 'Off';
end

% Under construction (not implemented yet)
if ~isfield(options, 'cores')
    % Default
    options.cores = 1;
end

%% Initialize

% Look for the folder in the data directory
options.dirFolder = dir(options.dataDir);    % The folders within the directory

% Create folder in which the output files will be stored
[~,~,~] = mkdir(options.outputDir);

%% Go through folders: load, process and save

if strcmp(options.singleFolder,'On')
    tFolders = 1;
else
    tFolders = length(options.dirFolder);
end

% Go through every folder and look for compress the state files.
for iFolder = 1:tFolders
    
    % If the selected content is a folder, and not one of those freaking
    % "." or ".."
    if (options.dirFolder(iFolder).isdir == 1 &&...
            ~strcmp(options.dirFolder(iFolder).name,'.') &&...
            ~strcmp(options.dirFolder(iFolder).name,'..') &&...
            strcmp(options.singleFolder,'Off')) ||...
            strcmp(options.singleFolder,'On')

        % If only one folder is given (So not a folder with subfolder with
        % the files but just one folder with files inside)
        if strcmp(options.singleFolder,'On')
            dir_file = dir(options.dataDir);
        else
            % Files in the folder
            dir_file = dir([options.dataDir '/' options.dirFolder(iFolder).name]);
            % Progress report
            fprintf('folder: %s\n',options.dirFolder(iFolder).name);
        end
        
        % Clear the temperary structure
        clear('endState');
        
        % Get structure
        try
            if strcmp(options.data.name{1},'All')
                [endState,meanState] = create_structures_from_files_all_fields;
            else
                [endState,meanState] = create_structures_from_files;
            end
            
            % Sort structure in order to make it causal
            sort_structure;
            
        catch
            warningMessage = sprintf('No state file constructed for folder %s. It could be that compressed state file already exists.     ',[options.dataDir '/' options.dirFolder(iFolder).name]);
            warning(warningMessage);
        end
        
        % Steps for saving output
        try
            if strcmp(options.singleFolder,'Off')
                % Get the subfolder of data
                if strcmp(options.subFolder,'On')
                    subFolder = [options.dirFolder(iFolder).name '/'];
                    %             subFolder = strsplit(options.dirFolder(i).name, '/');
                    %             subFolder = [subFolder{end-1} '/'];
                else
                    subFolder = '';
                end
                % Make folder
                [~,~,~] = mkdir([options.outputDir '/' subFolder]);
                
                % Save
                save([options.outputDir '/' subFolder 'point_' options.outputName '_' options.dirFolder(iFolder).name], '-struct', 'endState')
                save([options.outputDir '/' subFolder 'mean_' options.outputName '_' options.dirFolder(iFolder).name], '-struct', 'meanState')
            elseif strcmp(options.singleFolder,'On')
                save([options.outputDir '/' 'point_' options.outputName], '-struct', 'endState')
                save([options.outputDir '/' 'mean_' options.outputName], '-struct', 'meanState')
            end
            
            
        catch
            fprintf('No file saved for: %s \n',options.dirFolder(iFolder).name);
        end
    end
end

%% Sub-functions

    function varargout = create_structures_from_files
        % Get the compressed states of the folder
        
        l = 1;  % frame counter
        
        % Check every file of folder
        for j = 3:length(dir_file)
            
            % If the folder does contain the tag, then process
            if strfind(dir_file(j).name,options.tag) > 0
                
                % load the file
                state = load([options.dataDir '/' options.dirFolder(iFolder).name '/' dir_file(j).name]);
                
                % Check for all the different desired fields
                for k = 1:length(options.data.name)
                    
                    % Try to assume that the field exists
                    try
                        
                        % Try to process it  in the way most of the fields
                        % are constrcuted
                        try
                            
                            % Get the last time point of the file
                            endState.(options.data.name{k}).counts(:,:,l) = state.(options.data.name{k}).counts(:,:,end);
                            meanState.(options.data.name{k}).counts(:,:,l) = mean(state.(options.data.name{k}).counts(:,:,:),3);
                        catch
                            
                            % The execptions on the rule need a different
                            % approach
                            switch options.data.name{k}
                                case 'Time'
                                    endState.Time.values(1,1,l) = state.Time.values(:,1,end);
                                case 'Mass'
                                    endState.Mass.total(1,:,l) = state.Mass.total(1,:,end);
                                otherwise
                                    fprintf('no data saved for: %s',options.data.name{k});
                            end
                            
                        end
                        
                    catch
                        fprintf('no data saved for: %s\n',dir_file(j).name);
                    end
                end

                % New time frame
                l = l+1;
            end
            
        end
        varargout{1} = endState;
        varargout{2} = meanState;
    end

    function varargout = create_structures_from_files_all_fields
        % Get the all the compressed states of the folder
        
        l = 1;  % frame counter
        
        % Check every file of folder
        for j = 3:length(dir_file)
            
            % If the file does contain the tag, then process
            if strfind(dir_file(j).name,options.tag) > 0
                
                % load the file
                state = load([options.dataDir '/' options.dirFolder(iFolder).name '/' dir_file(j).name]);
                
                % Get the field names
                fieldNames = fieldnames(state);
                
                % Get all the data every field
                for k = 1:length(fieldNames)
                    subFieldNames = fieldnames(state.(fieldNames{k}));
                    for m = 1:length(subFieldNames)
                        try
                            % Get the last value of every data type
                            endState.(fieldNames{k}).(subFieldNames{m})(:,:,l) = state.(fieldNames{k}).(subFieldNames{m})(:,:,end);
                            meanState.(fieldNames{k}).(subFieldNames{m})(:,:,l) = mean(state.(fieldNames{k}).(subFieldNames{m})(:,:,:),3);
                        catch
                            
                        end
                        
                    end
                    
                end
                
                l = l + 1;
            end
            
        end
        
        varargout{1} = endState;
        varargout{2} = meanState;
        
    end

    function sort_structure
        % Sort the values in order to have a causal timeline
        
        % Get the right order
        [~,order] = sort(endState.Time.values(:));
        
        % Order all the  fields
        for k = 1:length(options.data.name)
            fieldNames = fieldnames(endState);
            for n = 1:length(fieldNames)
                subFieldNames = fieldnames(endState.(fieldNames{n}));
                for m = 1:length(subFieldNames)
                    try
                        endState.(fieldNames{n}).(subFieldNames{m}) = endState.(fieldNames{n}).(subFieldNames{m})(:,:,order);
                        meanState.(fieldNames{n}).(subFieldNames{m}) = meanState.(fieldNames{n}).(subFieldNames{m})(:,:,order);
                    catch
                        
                    end
                    
                end
            end
        end
        
    end
end
