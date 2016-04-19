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
%   Example:
%   structure = find_files('path/to/folder','FileNamePart')

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 06/01/2016

function out = find_files(varargin)

% Inputs
directory = varargin(1);
tag = varargin{2};

if nargin == 3
    showProgress = varargin{3};
else
    showProgress = 'Off';
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
            i = i + 1;
        end
    end

    function lookfortagedfiles
        % Search every folder to find the folders with the right tag
        k = 1;
        for i = 1:length(folderList)
            files = dir(folderList{i});
            if strcmp(showProgress,'On')
               display_progress(i,length(folderList)); 
            end
            for j = 3:length(files)
                if strfind(files(j).name,tag) > 0
                    % Save file name
                    out(k).file = files(j).name;
                    % Save directory name
                    out(k).folder = folderList{i};
                    k = k + 1;
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
