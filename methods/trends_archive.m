%trends_archive Analyse the trends of the parameters of the summary
%   file.
%   archive = trends_archive(archive).
%   The selected parameters are fitted using a first order polynomial fit
%   (polyfit(x,y,1)).
%
%   Generated fields:
%   archive.set.trends
%   archive.set.simulation.trends
%
%   Required fields:
%   archive.set.simulation.folder
%
%   Required files:
%   summary.mat (Default)
%
%   'set'           - Only check one or multiple specific sets.
%                     Example: 'set',[1,2,10],...
%   'fileTypeTag'   - Select the tag that the folder should contain. The
%                     default is value is "summary" such that the summary
%                     files are analyzed.
%   'fields'        - Select the fields that should be fitted. Default:
%                     fields = {'proteins','rnas'}. Note: This function
%                     might not support a certain field. This function
%                     should be adjusted accordingly if this is the case.
%   'progressState' - Show progress of this function ('On' or 'Off').
%                     Default value is 'Off'.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/05/2016

function archive = trends_archive(varargin)

% Mandatory input arguments
archive = varargin{1};

% -- Options and adjustments --%

% Initiate structure
if nargin > 1;
    options = struct(varargin{2:end});
else
    options = struct;
end
if ~isfield(options,'fileTypeTag')
    % Find file with name:
    options.fileTypeTag = 'summary.mat';
end
if ~isfield(options,'fields')
    % Analyse fields:
    options.fields = {'proteins','rnas','mass'};
end
if ~isfield(options,'progressState')
    % Progress update 'On' or 'Off':
    options.progressState = 'On';
end
if ~isfield(options,'set')
    % Total number of strains
    tStrain = length(archive.set);
    options.set = 1:tStrain;
else
    tStrain = length(options.set);
end
if ~isfield(options,'redo')
    options.redo = 'Off';
end
if ~isfield(options,'time_target')
    options.time_target = [0,1,2,3,4,5,6,7]*10000;
end

% -- Execution -- %

% Get whether the fileTypeTag is direct or not
C = strsplit(options.fileTypeTag,'.');
if strcmp(C(end),'mat')
   options.direct = 'On'; 
else
    options.direct = 'Off'; 
end

% Get updates on new line
if strcmp(options.progressState,'On')
    fprintf('Progress:\n     ')
end

% For every strain
for iStrain = options.set
    
    % Display progress
    if strcmp(options.progressState,'On')
        display_progress(iStrain,tStrain)
    end
    
    % Add field if not exist yet
    if ~isfield(archive.set(iStrain).simulation(1),'trends')
        archive.set(iStrain).simulation(1).trends = [];
    end
    
    % For every simulation
    for iSimulation = 1:length(archive.set(iStrain).simulation)
        
        % Only process when it has not been processed before or everything
        % should be redone.
        if isempty(archive.set(iStrain).simulation(iSimulation).trends) || strcmp(options.redo,'On')
            
            % Clear temperary variable
            clear temp
            
            % Clear trends
            archive.set(iStrain).simulation(iSimulation).trends = [];
            
            % Get the path to the simulation file
            if strcmp(options.direct,'On')
                folder = archive.set(iStrain).simulation(iSimulation).folder;
                summaryPath = {[folder '/' options.fileTypeTag]};
            else
                summaryPath = get_files(archive,options.fileTypeTag,iStrain,iSimulation);
            end
            
            % Load summary file (only the required fields)
            summary = load(summaryPath{1},options.fields{:},'time');
            
            % -- Add quantaties -- %
            
            % Time
            time = summary.time;
            
            % Proteins
            if sum(strcmp(options.fields,'proteins')) == 1
                % Total proteins
                data = double(summary.proteins(1,:));
                % Fit the data
                P = polyfit(time,data,1);
                % Add values
                archive.set(iStrain).simulation(iSimulation).trends_fit.proteins = P;
                % initiate previous time point
                previousTimePoint = 0;
                % Intermediate time points
                for iTimeTarget = 1:length(options.time_target)
                    [value,timePoint]=min(abs(time-options.time_target(iTimeTarget)));
                    if timePoint == previousTimePoint
                        % Don't reproduce value
                       break 
                    end
                    archive.set(iStrain).simulation(iSimulation).trends(iTimeTarget).time = timePoint;
                    archive.set(iStrain).simulation(iSimulation).trends(iTimeTarget).proteins = data(timePoint);
                    if value > 1
                        break
                    end
                    previousTimePoint = timePoint;
                end
            end
            
            % Rnas
            if sum(strcmp(options.fields,'rnas')) == 1
                % Total rnas
                data = double(summary.rnas(1,:));
                % Fit the data
                P = polyfit(time,data,1);
                % Add values
                archive.set(iStrain).simulation(iSimulation).trends_fit.rnas = P;
                % initiate previous time point
                previousTimePoint = 0;
                % Intermediate time points
                for iTimeTarget = 1:length(options.time_target)
                    [value,timePoint]=min(abs(time-options.time_target(iTimeTarget)));
                    if timePoint == previousTimePoint
                        % Don't reproduce value
                       break 
                    end
                    archive.set(iStrain).simulation(iSimulation).trends(iTimeTarget).time = timePoint;
                    archive.set(iStrain).simulation(iSimulation).trends(iTimeTarget).rnas = data(timePoint);
                    if value > 1
                        break
                    end
                    previousTimePoint = timePoint;
                end
            end
            
            % Mass
            if sum(strcmp(options.fields,'mass')) == 1
                % Total rnas
                data = double(summary.mass(1,:));
                % Fit the data
                P = polyfit(time,data,1);
                % Add values
                archive.set(iStrain).simulation(iSimulation).trends_fit.mass = P;
                % initiate previous time point
                previousTimePoint = 0;
                % Intermediate time points
                for iTimeTarget = 1:length(options.time_target)
                    [value,timePoint]=min(abs(time-options.time_target(iTimeTarget)));
                    if timePoint == previousTimePoint
                        % Don't reproduce value
                       break 
                    end
                    archive.set(iStrain).simulation(iSimulation).trends(iTimeTarget).time = timePoint;
                    archive.set(iStrain).simulation(iSimulation).trends(iTimeTarget).mass = data(timePoint);
                    if value > 1
                        break
                    end
                    previousTimePoint = timePoint;
                end
            end
            
            % Add new cases below if needed.
            
        end
        
    end
    
    %-- Get an average value for the fields --%
    
    % Field names
    nameFields = fields(archive.set(iStrain).simulation(1).trends);
    
    % For every trend
    for iField = 1:length(nameFields)
        
        % Total number of simulations
        tSim = length(archive.set(iStrain).simulation);
        
        % clear values
        values = [];
        
        % Determine values
        for iSim = 1:tSim
            values(:,iSim)=[archive.set(iStrain).simulation(iSimulation).trends.(nameFields{iField})];
        end
        
        % Calculate avarage
        average = nanmean(values,2);
        
        % Number of trends
        tTrend = length(average);
        
        % Add value to archive
        for iTrend = 1:tTrend
            archive.set(iStrain).trends(iTrend).(['average_' nameFields{iField}]) = average(iTrend);
        end
        
    end
    
end

end