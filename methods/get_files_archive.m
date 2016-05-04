%Get the paths to all the documents of a certain strain of a certain type. 
%
%   fileDir = get_files(archive,tag,strainNumber)
%
%   fileDir: Character array with all the paths to the files
%   archive: The archive 
%   tag: A part of the name of the file you are looking for, which should
%   be present in the folders that are given in the archive of that specific
%   strain.
%   strainNumber: The number of the strain according to archive.set.
%   
%   Optional:
%   'warning'   -'On' or 'Off', supress warning statements of this function.
%                 Default: 'On'
%   'multi'     -'On' or 'Off', multiple files are allowed in a folder and
%                 will be added to the archive. Default: 'Off'

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/04/2016

function varargout = get_files_archive(varargin)

% Mandatory input parameters
archive = varargin{1};
fileTypeTag = varargin{2};
strain = varargin{3};

% Optional input paramters
warningState = 'On';
multiState = 'On';
simulationNumber = 0; % All simulations

if nargin > 3;
    for i = 4:nargin
        if strcmp(varargin{i},'warning')
            warningState = varargin{i+1};
        end
        if strcmp(varargin{i},'multi')
            multiState = varargin{i+1};
        end
        if strcmp(varargin{i},'simulation')
            simulationNumber = varargin{i+1};
        end
    end
end

% Preallocation
fileDir = {};

% Get all files of all simulations if the simulation number is not
% specified
if simulationNumber == 0 
    simulationNumber = 1:length(archive.set(strain).simulation);
end

for i = simulationNumber
       
        % Find file that contains the tag
        tagFile = find_files(archive.set(strain).simulation(i).folder,fileTypeTag);
        % If nothing is found, change the tagFile output
        if length(tagFile) > 1 && strcmp(multiState,'Off')
            % Only one file can be used if multi mode is 'Off'
            errorTest = sprintf('Too many files in the folder "%s" matches the tag "%s"',...
                archive.set(strain).simulation(i).folder,fileTypeTag);
            error(errorTest)
        elseif length(tagFile) > 1 && strcmp(multiState,'On')
                fileDir(end+1,1) = {[archive.set(strain).simulation(i).folder '/' tagFile(1).file]};
            for fileNumber = 2:length(tagFile)
                % Safe path
                fileDir(end,fileNumber) = {[archive.set(strain).simulation(i).folder '/' tagFile(fileNumber).file]};
            end
        elseif isempty(tagFile)
            if strcmp(warningState,'On')
                % Mention when a file is missing
                warningMessage = sprintf('No "%s" file found in folder "%s"\n',fileTypeTag,archive.set(strain).simulation(i).folder);
                warning(warningMessage)
                fileDir(end+1,1) = {nan};
            end
        elseif length(tagFile) == 1
            % As it should be (just one file that is selected that contains
            % the information about the simulation)
            fileDir(end+1,1) = {[archive.set(strain).simulation(i).folder '/' tagFile(1).file]};
        end
%     end
end

% Get output arguments
varargout{1} = fileDir;

end
