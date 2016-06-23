%Plot a quantity of all or a set of simulations of a set.
%
%   In case of using data from state files:
%   plot_archive(archive,strainNumber,{field,subfield},molecule,location)
%   Note: molecule and location must still be defined when plotting other
%   quantities such as mass, geometry and etc. In that case, molecule and
%   location are typically "1". This corresponds to the location according
%   to the state files. 
%   
%   In case of using data from summary files:
%   plot_archive(archive,strainNumber,{field},location)
%   Note: location should be "1" if there is just 1 location
%
%   archive: is the archive. 
%   strainNumber: is the entry number of the targeted strain. 
%   field: The type of quantity to plot
%   subfield: The specific quantity of the field.
%   molecule: Typically the number of the molecule to plot. Might be
%   different for quantities such as "Mass" and "Geometry".
%   location: Typically the location of the molecule. The number are 1:6
%   where number "1" is the cytoplasma.
%
%   Created fields:
%   None
%
%   Required fields:
%   archive.set.simulation.folder
%
%   Required files:
%   data files like: state file or summary file
%   
%   Examples: 
%   plot_archive(archive,9,{'Metabolite','counts'},50,1)
%   plot_archive(archive,9,{'mass'},1)
%
%   Optional:
%   'figureNumber'  - Set the figure number of the plot. Default generates a
%                     new figure.
%   'fileTypeTag'   - Target file that should be used for this analysis.
%                     This function looks for files that contain a certain
%                     word (or string) in the file name.
%   'simulation'    - Select simulations to plot (example [1,2,5])
%   'symbol'        - Set the symbol for the data points (example 'o')
%   'step'          - Set the step size. Meaning, how many data points
%                     should be skipped when plotting (example 10)
%   'glueState'     - Glue togheter the multiple state files when needed.
%                     'glueState' can either be 'mean' or 'end'. meaning
%                     the mean or last value of the state file.
%   'progress'      - Show progress. ('On' or 'Off')
%   'normalize'     - Normalize values with max value

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/04/2016

function plot_archive(varargin)


%% Mendatory input arguments
archive = varargin{1};
strainNumber = varargin{2};
if length(varargin{3}) == 1
    fileType = 'summary';
    field = varargin{3}{1};
    location = varargin{4};
    fileTypeTag = 'summary.mat';
elseif length(varargin{3}) == 2
    fileType = 'state';
    field = varargin{3}{1};
    subfield = varargin{3}{2};
    molecule = varargin{4};
    location = varargin{5};
    fileTypeTag = 'mean_compressed_state.mat';
end

plot_options = {};

%% Adjust default settings

% Set figure number
entry = find(strcmp(varargin,'figureNumber'));
if ~isempty(entry)
    figureNumber = varargin{entry+1};
    figure(figureNumber);
else
    figure();
end
% Targeted files for information
entry = find(strcmp(varargin,'fileTypeTag'));
if ~isempty(entry)
    fileTypeTag = varargin{entry+1};
else
end
% Select specific simualtions
entry = find(strcmp(varargin,'simulation'));
if ~isempty(entry)
    simulationTargets = varargin{entry+1};
else
end
% Symbol
entry = find(strcmpi(varargin,'symbol'));
if ~isempty(entry)
    plot_options(end+1) = {varargin{entry+1}};
else
    plot_options(end+1) = {'-'};
end
% Color
entry = find(strcmpi(varargin,'color'));
if ~isempty(entry)
     plot_options(end+1) = {'Color'};
     plot_options(end+1) = {varargin{entry+1}};
end
% Select step size
entry = find(strcmp(varargin,'step'));
if ~isempty(entry)
    step = varargin{entry+1};
else
    step = 1;
end
% Select glueState
entry = find(strcmp(varargin,'glueState'));
if ~isempty(entry)
    glueState = varargin{entry+1};
else
    glueState = 'end';
end
% Progress state
entry = find(strcmp(varargin,'progress'));
if ~isempty(entry)
    options.progress = varargin{entry+1};
else
    options.progress = 'On';
end
% Normalize the values
entry = find(strcmp(varargin,'normalize'));
if ~isempty(entry)
    options.normalize = varargin{entry+1};
else
    options.normalize = false;
end

%% Processing

% Get whether the fileTypeTag is direct or not
C = strsplit(fileTypeTag,'.');
if strcmp(C(end),'mat')
   options.direct = 'On'; 
else
    options.direct = 'Off'; 
end

% Total number of simulations
tSimulations = length(archive.set(strainNumber).simulation);

% Find the desired files
fileDir = {};
if strcmp(options.direct,'Off')
    fileDir = get_files(archive,fileTypeTag,strainNumber,'multi','On');
elseif strcmp(options.direct,'On')
    for iSim = 1:tSimulations
        folder = archive.set(strainNumber).simulation(iSim).folder;
        fileDir(end+1,1) = {[folder '/' fileTypeTag]};
    end
end

% Maximal number of files associated with one simulation
tSimFiles = size(fileDir,2);
if tSimFiles > 1
    multiFiles = 'On';
    fprintf(['In order to speed up the plotting of data, it''s adviced to \n'...
        'make a compressed file of the state files and use these files to  \n'...
        'plot. These files can be constructed with xxx and selected by setting  \n',...
        '''fileTypeTag'' to ''point_compressed'' or ''mean_compressed  \n'''...
        'accordingly.\n'...
        'Example: archive.plot(1,{''Metabolite'',''counts''},50,1,''fileTypeTag'',''point_compressed'')\n']);
else 
    multiFiles = 'Off';
end

% Select simulations
if ~exist('simulationTargets','var')
    simulationTargets = 1:tSimulations;
end

% Progress spacer
if strcmp(options.progress,'On')
    fprintf('Progress:\n     ');
end

% Load and plot data
for iSimulation = simulationTargets
    
    % Display progress
    if strcmp(options.progress,'On')
        display_progress(iSimulation,length(simulationTargets))
    end
    
    try
        
    % Load database
    if strcmp(fileType,'state')
        % Check if multiple files should be glued together
        if strcmp(multiFiles,'Off')
            % For single state files
            state = load(fileDir{iSimulation,1},field,'Time');
        elseif strcmp(multiFiles,'On')
            % For when there are multiple state files that should be glued
            % together
            % For every state file
            for iSimFile = 1:tSimFiles
                if ~isempty(fileDir{iSimulation,iSimFile})
                    % If state file exist, load file
                    tempState = load(fileDir{iSimulation,iSimFile},field,'Time');
                    if strcmp(glueState,'mean')
                        % Get the mean values of the state file
                        state.(field).(subfield)(:,:,iSimFile)=mean(tempState.(field).(subfield),3);
                        state.Time.values(:,:,iSimFile) = round(mean(tempState.('Time').('values'),3));
                    elseif strcmp(glueState,'end')
                        % Get the max values of the state file
                        state.(field).(subfield)(:,:,iSimFile)=max(tempState.(field).(subfield),[],3);
                        state.Time.values(:,:,iSimFile) = max(tempState.('Time').('values'),[],3);
                    end
                end
            end
            % sort the values from early to late in the simulation
            [~,order] = sort(state.Time.values);
            state.(field).(subfield) = state.(field).(subfield)(:,:,order);
            state.Time.values = state.Time.values(:,:,order);
        end
    elseif strcmp(fileType,'summary')
        % For summary files
        state = load(fileDir{iSimulation,1},field,'time');
    end
    
    % Get quantities
    if strcmp(fileType,'state')
        % For state files
        time = squeeze(state.Time.values);
        data = state.(field).(subfield)(molecule,location,:);
        data = shiftdim(data,5)';
    elseif strcmp(fileType,'summary')
        % For summary files
        time = squeeze(state.time);
        data = squeeze(state.(field)(location,:));
    end
    
    % Normalize data
    if options.normalize
       data = data/max(max(data)); 
    end
    
    % Plot data
    hold on
    plot(time(1:step:end),data(:,1:step:end),plot_options{:});
    hold off
    
    catch
        warningMessage = sprintf('Cannot plot the data fields %s:%s of simulation %d set %d',...
            field,subfield,iSimulation,strainNumber);
       warning(warningMessage) 
    end
end

% Layout
if strcmp(fileType,'state')
    % For state files
    title([field ': ' subfield])
elseif strcmp(fileType,'summary')
    % For summary files
    title([field])
end

% Set x label
xlabel('Time')

end
