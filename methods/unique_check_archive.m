%unique_archive_list    Checks the uniqueness of the simulations of every 
%   strain by comparing their seeds.
%   list = viable_archive_list(list)
%
%   Generated fields: 
%   list.simulation.unique
%   
%   Required fields: 
%   list.simulation.seed
%
%   Required files: 
%   None
%
%   Optional:
%   'remove'    - Set whether the non unique simulations should be removed
%                 or not. The default setting is 'On'.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 01/26/2016

function archive = unique_check_archive(varargin)

% get input
archive = varargin{1};

% Default setting
removeState = 'On'; % remove non unique simulations from list

% Adjust default setting
for i = 2:nargin
    if strcmp(varargin{i},'remove')
        removeState = varargin{i+1};
    end
end

non_unique_counter = 0;

% Check every strain
for i = 1:length(archive.set)
    
    % Check every simulation file
    for j = 1:length(archive.set(i).simulation)
        
        % Check if seed is unique
        if isempty(find([archive.set(i).simulation(1:j-1).seed] == archive.set(i).simulation(j).seed))
            % If unique
            archive.set(i).simulation(j).unique = 1;
        else
            % If not unique
            archive.set(i).simulation(j).unique = 0;
            non_unique_counter = non_unique_counter+1;
        end
    end
end

% Remove non unique states
if strcmp(removeState,'On')
    % Check every strain
    for i = 1:length(archive.set)
        % Keep only the unique simulations
        archive.set(i).simulation = archive.set(i).simulation([archive.set(i).simulation.unique]==1);
    end
end

fprintf('%d non unique simulation(s) found\n',non_unique_counter)

end
