%Find a file with a certain tag in its name. Files in both the folder 
%   and its subfolders are checked. the output is a structure containing 
%   the path and the name of the folder. The fields of the output structure 
%   are file and folder. Both inputs "directory" and "tag" should be a 
%   string.
%
%   out = find_files(directory,tag)
%   
%   Where directory and tag are both string values and out is a structure
%   with the fields "file" and "folder".
%
%   'short'     - Setting 'short' to save con significantly speed up the
%                 find folder process. This option should only be used when
%                 all the simulation files are at the same level!!! This
%                 because the find folder project will stop when a folder
%                 is found without any simulations inside.
%
%   Examples:
%   structure = find_files('path/to/folder','FileNamePart')
%   out = find_files(directory,tag,'Off') % show no progress
%   out = find_files(directory,tag,'short')

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 06/01/2016

function out = find_files(varargin)

% Inputs
directory = varargin(1);
tag = varargin{2};

if nargin > 2
    if strcmpi(varargin{3},'On') || strcmpi(varargin{3},'Off')
        showProgress = varargin{3};
    else
        showProgress = 'Off';
    end
else
    showProgress = 'Off';
end

% shortState: 'On'. Stop looking for folders at the point when no folders
% are being found.
iEntry = find(strcmp(varargin(3:end),'short'));
if ~isempty(iEntry)
    shortState = 'On';
else
    shortState = 'Off';
end

% Get whether the fileTypeTag is direct or not
C = strsplit(tag,'.');
if strcmp(C(end),'mat')
   directState = 'On'; 
else
   directState = 'Off'; 
end

% Global variable
folderList = directory;

if strcmp(showProgress,'On')
    fprintf('Progress get folders:      ');
end

% Find all subfolders
lookforallfolders

if strcmp(showProgress,'On')
    fprintf('Progress get files:      ');
end

% Find all files with the tag in all folders
lookfortagedfiles

    function lookforallfolders
        % Make a list of all subfolders in the directory
        i = 1;
        while i <= length(folderList)
            
            if strcmp(showProgress,'On')
               display_progress(i,length(folderList)); 
            end
            
            % Find all files
            files = dir(folderList{i});
            % Check if they are folders
            dirFlags = [files.isdir];
            % Only select the folders
            subFolders = files(dirFlags);
            for j = 3:length(subFolders)
                % Add all new folders to folder list
                folderList{end+1}=[folderList{i} '/' subFolders(j).name];
            end
            % Stop if only '.' or '..' is found.
            if (strcmp(shortState,'On') && length(subFolders) < 3)
                display_progress(1,1); 
                break
            end
            
            i = i + 1;
        end
    end

    function lookfortagedfiles
        % Search every folder to find the folders with the right tag
        
        out = [];
        tFolder = length(folderList);
        
        % If the target is specific (example: options.mat)
        if strcmpi(directState,'On')
            for i = 1:tFolder
                if strcmp(showProgress,'On')
                    display_progress(i,tFolder);
                end
                state = exist([folderList{i} '/' tag],'file');
                if state > 0
                    % Save file name
                    out(end+1).file = tag;
                    % Save directory name
                    out(end).folder = folderList{i};
                end
            end
        end
        
        % if the tag is not specific (example: options)
        if strcmpi(directState,'Off')
            for i = 1:tFolder
                files = dir(folderList{i});
                for j = 3:length(files)
                    if strfind(files(j).name,tag) > 0
                        % Save file name
                        out(end+1).file = files(j).name;
                        % Save directory name
                        out(end).folder = folderList{i};
                    end
                end
            end
        end
    end

% Output an empty output if nothing is found
if ~exist('out','var')
    out = struct([]);
    warningMessage = sprintf('No files found in the folder %s with the tag %s',directory{1},tag);
    warning(warningMessage);
end

end
