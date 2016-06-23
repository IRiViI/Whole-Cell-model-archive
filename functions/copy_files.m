%Copy target files to another location. This function is designed in order
%   to relocate files of the whole cell model.
%   copy_files_archive(path_data,path_destination)
%   Copy all the default files within subfolders of the folder path_data to 
%   the new location path_destination. The folder structure in preserved.
%   
%   Example:
%   copy_files('path/to/folder','path/to/destination',...
%       'documentTypes',{'summary.mat','options.mat'});
%   Move the files 'summary.mat' and 'options.mat' to their new location.
%   The subfolder structure is preserved. 

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/13/2016

function copy_files(varargin)

% Mandatory inputs
options.moveFromDir = varargin{1};
options.moveToDir = varargin{2};

% Default options
options.showProgress = 'On';
options.shortState = 'Off';
options.documentTypes =  {...
    'summary.mat',...
    'options.mat',...
    'parameters.mat',...
    'point_compressed_state.mat',...
    'mean_compressed_state.mat'};

% Adjust options
inputOptions = struct(varargin{3:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

% Process options
tDocument = length(options.documentTypes);

% Obtain folder structure
[~,lSubFolder] = lookforallfolders(options.moveFromDir,...
    'showProgress',options.showProgress,'shortState',options.shortState);

% Total number of sub folders
tSubFolder = length(lSubFolder);

% Progress
if strcmpi(options.showProgress,'On')
    fprintf('Progress:/n      ')
end

% For every found sub folder
for iSubFolder = 1:tSubFolder
    
    % Progress 
    if strcmpi(options.showProgress,'On')
        display_progress(iSubFolder,tSubFolder);
    end
    
    % Folder to check
    fromFolder = [options.moveFromDir '/' lSubFolder{iSubFolder}];
    
    % Folder to move file to
    toFolder = [options.moveToDir '/' lSubFolder{iSubFolder}];
    
    % Check every document type in folder
    for iDocument = 1:tDocument
        
        % File to check
        file = [fromFolder '/' options.documentTypes{iDocument}];
        
        % Check if file exist
        if exist(file,'file')
            
            % Make folder for file (if required)
            [~,~,~] = mkdir(toFolder);
            
            % Copy file
            copyfile(file ,toFolder);
            
        end
        
    end
    
end

end