%Determine the growth rate of cells and strains. Additionally is the mass
%   at the end of the simulation saved.
%
%   archive = growth_rate_archive(archive)
%
%   Generated fields:
%   archive.set.simulation.growthRate.(fit,std,mean,mass.end)
%   archive.set.growthRate.(averageFit,averageStd,averageMean)
%   archive.set.simulation.mass
%
%   Required fields:
%   archive.set.simulation.folder
%
%   Required files:
%   Compressed state files
%
%   Optional:
%
%   'set'           - Only check one or multiple specific sets.
%                     Example: 'set',[1,2,10],...
%   'progress'      - Display the progress of this function. Default
%                     setting is 'Off'
%   'plot'          - Plot the actual data and the fit while running this
%                     function. Default setting is 'Off'.
%   'fileTypeTag'   - Target file that should be used for this analysis.
%                     This function looks for files that contain a certain
%                     word (or string) in the file name.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/05/2016

function archive = growth_rate_archive(varargin)

% Mendatory inputs parameter
archive = varargin{1};

% Default settings
progressState = 'On';
plotState = 'Off';
fileTypeTag = 'point_compressed_state.mat';
location = [1,2,4,5,6];
sets = 1:length(archive.set);
options.redo = 'Off';

% Change settings
if nargin > 1;
    for i = 2:nargin
        if strcmp(varargin{i},'progress')
            progressState = varargin{i+1};
        end
        if strcmp(varargin{i},'plot')
            plotState = varargin{i+1};
        end
        if strcmp(varargin{i},'fileTypeTag')
            fileTypeTag = varargin{i+1};
        end
        if strcmp(varargin{i},'set')
            sets = varargin{i+1};
        end
        if strcmp(varargin{i},'redo')
            options.redo = varargin{i+1};
        end
    end
end

% Add settings of the function
archive.settings.growthRate.location = location;
archive.settings.growthRate.output = {'P(1)*t+P(2)'};

% Progress spacer
% Get updates on new line
if strcmp(progressState,'On')
    fprintf('Progress:\n     ')
end

% Total number of strains
tStrain = length(sets);

% For every strain
for strainNumber = sets
    
    % Progress
    if strcmp(progressState,'On')
        display_progress(strainNumber,tStrain);
    end
    
    % Get the location of the desired files
    %     fileDir = get_files(archive,fileTypeTag,strainNumber);
    
    % Number of simulations in set
    tSimulation = length(archive.set(strainNumber).simulation);
    
    % clean and preallocate temperary variables
    tempFit = zeros(tSimulation,2);
    tempStd = zeros(tSimulation,1);
    tempMean = zeros(tSimulation,1);
    
    % Add field if not exist yet
    if ~isfield(archive.set(strainNumber).simulation(1),'growth_rate')
        archive.set(strainNumber).simulation(1).growth_rate = [];
    end
    
    % For very simultion
    for simulationNumber = 1:tSimulation
        
        % Only process when it has not been processed before or everything
        % should be redone.
        if isempty(archive.set(strainNumber).simulation(simulationNumber).growth_rate) || strcmp(options.redo,'On')
            
            try
            % Load the state file
            state = load_state(archive,strainNumber,simulationNumber,...
                'fields',{'Time','Mass'},...
                'fileTypeTag',fileTypeTag);
            
            % Get the values for the mass
            y = squeeze(sum(state.Mass.total(1,location,:),2));
            
            % Get time values
            t = squeeze(state.Time.values);
            
            % Determine the change of mass
            dy = y(2:end)-y(1:end-1);
            dt = t(2:end)-t(1:end-1);
            dydt = dy./dt;
            
            % Fit change of mass
            fit = polyfit(t(1:end-1),dydt,1);
            
            catch
                warningMessage = sprintf('No file found for set %d simulation %d\n',strainNumber,simulationNumber);
                warning(warningMessage)
                fit = [nan,nan];
                dydt = nan;
                y = nan;
            end
            % Safe values
            archive.set(strainNumber).simulation(simulationNumber).growthRate.fit = fit;
            tempFit(simulationNumber,:) = fit;
            archive.set(strainNumber).simulation(simulationNumber).growthRate.std = std(dydt);
            tempStd(simulationNumber) = std(dydt);
            archive.set(strainNumber).simulation(simulationNumber).growthRate.mean = mean(dydt);
            tempMean(simulationNumber) = mean(dydt);
            
            % Other parameter
            archive.set(strainNumber).simulation(simulationNumber).mass.end = y(end);
            archive.set(strainNumber).simulation(simulationNumber).mass.begin = y(1);
            
            % Plot Data
            if strcmp(plotState,'On')
                figureNumber = 99; % Why not
                figure(figureNumber);
                clf(figureNumber);
                figure(figureNumber)
                hold on
                plot(fit(2)+fit(1)*t(1:end-1))
                plot(dydt)
                hold off
                titleMessage = sprintf('Strain: %d, simulation: %d',strainNumber,simulationNumber);
                title(titleMessage)
                xlabel('time')
                ylabel('Mass')
                pause(1)
            end
            
        end
        
    end
    
    % Determine the average value for the whole strain
    archive.set(strainNumber).growthRate.averageFit = mean(tempFit,1);
    archive.set(strainNumber).growthRate.averageStd = mean(tempStd,1);
    archive.set(strainNumber).growthRate.averageMean = mean(tempMean,1);
end


end
