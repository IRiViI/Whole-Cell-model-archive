function varargout = lookforallfolders(folderList,varargin)
% Find al folders in folder list and add them to folderlist. The list grows
% every itteration.
%
%   'showProgress'  - Show progress ('On' or 'Off')
%   'shortState'    - Stop looking for a folder after a folder without
%                     subfolders is found. This options speeds up the
%                     foldering finding process drastically, however it
%                     only results in a complet set of folders if the
%                     folder structure is conistent. ('On' or 'Off')
%                     default 'Off'.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 06/13/2016

% Default options
options.showProgress = 'On';
options.shortState = 'Off';

% Adjust options
inputOptions = struct(varargin{1:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

% Process inputs
if ischar(folderList)
    folderList = {folderList};
end

% Make a list of all subfolders in the directory
subFolderList = {};
iFolder = 1;
while iFolder <= length(folderList)
    
    if strcmp(options.showProgress,'On')
        display_progress(iFolder,length(folderList));
    end
    
    % Find all files
    files = dir(folderList{iFolder});
    % Check if they are folders
    dirFlags = [files.isdir];
    % Only select the folders
    subFolders = files(dirFlags);
    for j = 1:length(subFolders)
        if ~strcmp(subFolders(j).name(1),'.')
            % Add all new folders to folder list
            folderList{end+1} = [folderList{iFolder} '/' subFolders(j).name];
        end
    end
    
    % Stop if only '.' or '..' is found.
    if (strcmp(options.shortState,'On') && length(subFolders) < 3)
        display_progress(1,1);
        break
    end
    
    % Update counter
    iFolder = iFolder + 1;
    
end

%% Set outputs

% Set folderList as output
varargout{1} = folderList;

% Make list of subFolders
tSymbol = length(folderList{1}); % Number of symbols of seed folder
tFolder = length(folderList); % Total number of folders
for iFolder = 1:tFolder
    tmp = folderList{iFolder}(tSymbol+2:end); % subfolder name
    if ~isempty(tmp) % If it's not the original folder
        subFolderList{iFolder} = tmp; % Add subfolder name
    else
        subFolderList{iFolder} = '';
    end
end

% Set subfolder output
varargout{2} = subFolderList;

end
