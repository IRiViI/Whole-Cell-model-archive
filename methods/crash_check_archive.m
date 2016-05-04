%Check if the simulation crashed. A simulation is considered to be crashed
%   when the simulation stopped running before 50000(default) seconds while
%   there is no pinching obsereved.
%   archive = crash_check_archive(archive)
%
%   Created fields:
%   archive.set.simulation.crash, where crashed has the value 0 or 1
%   respectively to not crashing and crashing.
%
%   Required fields:
%   archive.set.simulation.folder
%
%   Required files:
%   Compressed state file.
%
%   'set'           - Only check one or multiple specific sets.
%                     Example: 'set',[1,2,10],...
%   'fileTypeTag'   - Target file that should be used for this analysis.
%                     This function looks for files that contain a certain
%                     word (or string) in the file name.
%   'remove'        - Remove crashed simulations. Default is 'Off'.
%   'progress'      - Display the progress of this function. Default
%                     setting is 'On'.
%   'lifeTime'      - Set the max lifetime of simulation

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 02/16/2016

function archive = crash_check_archive(varargin)

% Mandatory input
archive = varargin{1};

% Targeted file
fileTypeTag = 'summary.mat';
removeState = 'Off';
progressState = 'On';
lifeTime = 50000;
sets = 1:length(archive.set);
options.redo = 'Off';

% Adjust default settings
if nargin > 1;
    for i = 2:nargin
        if strcmp(varargin{i},'fileTypeTag')
            fileTypeTag = varargin{i+1};
        end
        if strcmp(varargin{i},'remove')
            removeState = varargin{i+1};
        end
        if strcmp(varargin{i},'progress')
            progressState = varargin{i+1};
        end
        if strcmp(varargin{i},'lifeTime')
            lifeTime = varargin{i+1};
        end
        if strcmp(varargin{i},'set')
            sets = varargin{i+1};
        end
        if strcmp(varargin{i},'redo')
            options.redo = varargin{i+1};
        end
    end
end

% Settings usded to determine whether the simulation crashed or not.
if ~isprop(archive,'settings')
    archive.addprop('settings');
end
archive.settings.crash.lifeTime = lifeTime;

% Initiate progress line
if strcmp(progressState,'On')
    fprintf('Progress:\n      ')
end

% Check every strain
for strainNumber = sets
    
    % Display progress strain
    if strcmp(progressState,'On')
        display_progress(strainNumber,length(sets))
    end
    
    % Create field if it does not exist yet
    if ~isfield(archive.set(strainNumber).simulation(1),'crash')
        archive.set(strainNumber).simulation(1).crash = [];
    end
    
    % Number of simulations
    tSim = length(archive.set(strainNumber).simulation);
    
    % Check every simulation
    for simulationNumber = 1:tSim
        
        % Only process when it has not been processed before or everything
        % should be redone.
        if isempty(archive.set(strainNumber).simulation(simulationNumber).crash) ||...
                strcmp(options.redo,'On')
            
            % Load the state file
            if ~strcmp(fileTypeTag,'summary.mat')
                % state files
                state = archive.load_file(strainNumber,simulationNumber,...
                    'fileTypeTag',fileTypeTag,...
                    'fields',{'Geometry'});
                
                % summary files
            elseif strcmp(fileTypeTag,'summary.mat')
                state = archive.load_file(strainNumber,simulationNumber,...
                    'fileTypeTag',fileTypeTag,...
                    'fields',{'pinchedDiameter'});
            end
            
            % Check if the the simulation meet the criteria of crashing
            try
                % State files
                if state.Geometry.pinchedDiameter(end) ~= 0 &&...
                        archive.set(strainNumber).simulation(simulationNumber).trends(end).time < archive.settings.crash.lifeTime
                    
                    archive.set(strainNumber).simulation(simulationNumber).crash = 1;
                else
                    archive.set(strainNumber).simulation(simulationNumber).crash = 0;
                end
                
                % Summary files
            catch
                if state.pinchedDiameter(end) ~= 0 &&...
                        archive.set(strainNumber).simulation(simulationNumber).trends(end).time < archive.settings.crash.lifeTime
                    
                    archive.set(strainNumber).simulation(simulationNumber).crash = 1;
                else
                    archive.set(strainNumber).simulation(simulationNumber).crash = 0;
                end
            end
            
        end
        
    end
    
end


% Remove non unique states
if strcmp(removeState,'On')
    % Check every strain
    for iStrain = sets
        % Keep only the unique simulations
        archive.set(iStrain).simulation = archive.set(iStrain).simulation([archive.set(iStrain).simulation.crash]==0);
    end
end

end
