%Load a state file of the archive
%   state = load_state(archive,strainNumber)
%   Load the first simulation of strain "strainNumber"
%
%   state = load_state(archive,strainNumber,simulationNumber)
%   Specify a specific simulation to load.
%
%   'fileTypeTag'   - Change the type of state file to load.
%   'fields'        - Load specific field(s). Example: {'a','b'}.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/05/2016

function state = load_state(varargin)

% Process mandatory input argument
archive = varargin{1};
strain = varargin{2};

% Default settings
number = 1;
fileTypeTag = 'point_compressed_state.mat';
fieldSelect = {};

% Adjust default settings
for i = 3:nargin
    % If the third input is a number, then "number" should be set
    % accordingly
    if i == 3 && isnumeric(varargin{3});
        number = varargin{3};
    end
    if strcmp(varargin{i},'fileTypeTag')
        fileTypeTag = varargin{i+1};
    end
    if strcmp(varargin{i},'fields')
        fieldSelect = varargin{i+1};
    end
end

% Get whether the fileTypeTag is direct or not
options.direct = check_direct(fileTypeTag);

% Get the path to the file. It is faster when the fileTypeTag is precise
if ~options.direct
    % Get the file path
    fileDirRef = get_files(archive,fileTypeTag,strain);
    target = fileDirRef{number};
else % Direct
    % Get the file path
    target = [archive.set(strain).simulation(number).folder '/' fileTypeTag];
end

% Load state file
if isempty(fieldSelect)
    state = load(target);
else
    state = load(target,fieldSelect{:});
end

end
