%Copy target files to another location. This function is designed in order
%   to relocate files of the whole cell model.
%   copy_files_archive(path_data,path_destination)
%   Copy all the default files within subfolders of the folder path_data to 
%   the new location path_destination. The folder structure in preserved.
%   
%   copy_files(path_data,path_destination,target_documents)
%   The default of target_documents are the documents: 'summary.mat',
%   'options.mat' and 'parameters.mat'. This can be changed by specifying
%   the target_documents
%
%   Example:
%   copy_files('path/to/folder','path/to/destination', {'summary.mat',...
%   'options.mat'})
%   Move the files 'summary.mat' and 'options.mat' to their new location.
%   The subfolder structure is preserved. 

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/06/2016

function copy_files(varargin)

% Process inputs
if nargin >= 2
    data_dir = varargin{1};         % Data origin
    copy_folder = varargin{2};      % Data destination
else
    error('Not enough input arguments');
end

if nargin == 3
    target_document = varargin{3};  % Targets of the folder
else
    % Default
    target_document = {'summary.mat','options.mat','parameters.mat'};
end

% Create folder
[~,~,~] = mkdir(copy_folder);

% Look for the folder in the data directory
dir_folders = dir(data_dir);    % The folders within the directory

% Look for the files with the name of the target_document and copy them to
% the directory folder if possible. Also look at the subfolders for
% target_document.
for i = 3:length(dir_folders)
    try
        % Create folder
        [~,~,~] = mkdir([copy_folder '/' dir_folders(i).name]);
        for l = 1:length(target_document)
            % The path to the file (if exist)
            path_input = [data_dir '/' dir_folders(i).name '/' target_document{l}];
            % path to output
            path_output = [copy_folder '/' dir_folders(i).name '/' target_document{l}];
            % Copy file (if possible)
            copyfile(path_input ,path_output);
        end
    catch
    end
end
end
