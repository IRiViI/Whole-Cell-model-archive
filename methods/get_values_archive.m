function matrix = get_values_archive(archive,set,simulation,target,...
    element,location,varargin)
% Extract values of all the single state files and combine them into one
%   single matrix.
%   matrix = get_values_archive(archive,set,simulation,target,element,...
%       location)
%
%   Example: matrix = get_values_archive(archive,1,1,...
%       'Metabolite.counts',[1:10],[1:6]);
%       Get the information the the first 10 metabolites at all locations 
%       of the first simulation of the first set.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/27/2016

% Default options
options.fileTypeTag = 'state-';
options.stateStart = 0;                 
options.stateStop = inf;

% Adjust options
inputOptions = struct(varargin{1:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

% Get field and subfield of target
field = strsplit_archive(target,'.');

% Folder simulation
folder = archive.folder(set,simulation);

% Initiate state file counter
number = options.stateStart;

% Initiate structure
state = struct([]);

while true
    
    % path to next state file
    state_part_path = [folder '/' options.fileTypeTag num2str(number) '.mat'];
    
    % Load field of state file, if exist
    if exist(state_part_path,'file')
        state_part = load(state_part_path,field{1});
    else
        break
    end
        
    % Get correct elements of structure
    state(number+1).values = state_part.(field{1}).(field{2})(element,location,:);
    
    % Update state file counter (note: zero should be processed in while loop!)
    number = number + 1;
    
    % Stop if max is reached
    if number > options.stateStop
       break 
    end
    
end

% Combine elements to a single matrix
matrix = cat(3,state.values);

end