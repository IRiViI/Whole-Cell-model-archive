% Template for making new function for the whole cell model archive class.
%   ARCHIVE = TEMPLATE_ARCHIVE(ARCHIVE) run simulation with default
%   settings.
%   Stripped version of "template_archive.m".

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function varargout = sim_set_quantaties_archive(archive, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.fileTypeTag = 'mean_compressed_state.mat';
options.field = '';         % Mandatory input
options.subfield = '';      % Mandatory input
options.minTime = 1;
options.maxTime = inf;
options.groupSize = 60;
options.extractionProgress = false;
options.save = true;        % options: true or false
options.outFolder = 'default';
options.outName = 'default';

% Change settings

tmpOptions = struct(varargin{1:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

%% Process settings

% Determine which sets of the archive should be analyzed
if options.set == 0
    tSet = archive.sets;    % Number of sets to analyze
    options.set = 1:tSet;          % List of sets to analyze
end

% Set default save state, only safe is number of sets is one
if strcmpi(options.save,'default')
    if length(options.set) == 1
        options.save = true;
    else
        options.save = false;
    end
elseif islogical(options.save)
else
    error('Incorrect save setting')
end

% Determine save location and name
if options.save
    
    % Set save location
    if strcmpi(options.outFolder,'default')
        if length(options.set) == 1
            % If there is only one set selected, safe the file in the parent of the
            % folder of the first simulation
            if options.simulation == 0
                folder = archive.folder(options.set,1);
            else
                folder = archive.folder(options.set,options.simulation(1));
            end
            options.outFolder = [folder '/' 'output'];
        else
            datetime('now')
            % Else (when default), save in output folder
            options.outFolder = [archive.settings.dir '/' 'output'];
        end
    end
    
    % Set save name
    if strcmpi(options.outName,'default')
        if length(options.set) == 1
            % Simple name when only 1 set is used
            options.outName = ['analysis' '_' ...
                options.field '_' ...
                options.subfield '_' ...
                options.fileTypeTag];
        else
            % Add time information when multiple sets are used
            tmp = datestr(datetime('now'));
            time = strrep(tmp, ':', ' ');
            options.outName = ['analysis' '_' ...
                time '_' ...
                options.field '_' ...
                options.subfield '_' ...
                '_' options.fileTypeTag];
        end
    end
    
    % Create output folder
    [~,~,~] = mkdir(options.outFolder);
end

%% Run code

% Display progress
if options.progress
    fprintf('Initialisation: ');
end

% Get elements which should be checked (The elements which are non zero in
% all the simulations)
[blElement,nzseOptions] = archive.non_zero_state_elements('set',options.set,...
    'simulation',options.simulation,...
    'field',options.field,...
    'subfield',options.subfield,...
    'compartment',1,...
    'maxTime',options.maxTime,...
    'fileTypeTag',options.fileTypeTag);

% Update max Time point
if nzseOptions.maxTime < options.maxTime
    options.maxTime = nzseOptions.maxTime;
    warning('maxTime setting to small for one or more simulations. maxTime adjusted accordingly.')
end

% List of elements to check
lElement = find(blElement);
options.element = lElement;

% Total number of elements
tElement = length(blElement);
tCheckElement = length(lElement);

% Make element groups (In order to speed up the function)
tGroup = ceil(tCheckElement/options.groupSize);
lGroup = vec2mat(lElement,options.groupSize);

% Initiate output structures
meanValue = zeros(tElement,options.maxTime-options.minTime+1);
stdValue = zeros(tElement,options.maxTime-options.minTime+1);
varValue = zeros(tElement,options.maxTime-options.minTime+1);

% Initiate halve output structures
halve(1).meanValue = zeros(tElement,options.maxTime-options.minTime+1);
halve(1).stdValue = zeros(tElement,options.maxTime-options.minTime+1);
halve(1).varValue = zeros(tElement,options.maxTime-options.minTime+1);
halve(2).meanValue = zeros(tElement,options.maxTime-options.minTime+1);
halve(2).stdValue = zeros(tElement,options.maxTime-options.minTime+1);
halve(2).varValue = zeros(tElement,options.maxTime-options.minTime+1);

% Display progress
if options.progress
    fprintf('Done\n');
    fprintf('Progress:\n     ');
    display_progress(0,1)
end

% Update non zero elements
for iGroup = 1:tGroup
    
    % Get current element numbers
    nElement = lGroup(iGroup,:);
    nElement = nElement(nElement~=0); % Remove zero values
    
    % Extract element values of state files
    matrix = archive.sim_set_matrix('set',options.set,...
        'simulation',options.simulation,...
        'field',options.field,...
        'subfield',options.subfield,...
        'maxTime',options.maxTime,...
        'element',nElement,...
        'compartment',1,...
        'progress',options.extractionProgress);
    
    % Determine mean, std and var values
    meanValue(nElement,:) = nanmean(matrix,4);
    stdValue(nElement,:) = nanstd(matrix,[],4);
    varValue(nElement,:) = nanvar(matrix,[],4);
    
    % Determine of halve
    first = 1:ceil(size(matrix,4)/2);
    second = ceil(size(matrix,4)/2) + 1 : size(matrix,4);
    halve(1).meanValue(nElement,:) = nanmean(matrix(:,:,:,first),4);
    halve(1).stdValue(nElement,:) = nanstd(matrix(:,:,:,first),[],4);
    halve(1).varValue(nElement,:) = nanvar(matrix(:,:,:,first),[],4);
    halve(2).meanValue(nElement,:) = nanmean(matrix(:,:,:,second),4);
    halve(2).stdValue(nElement,:) = nanstd(matrix(:,:,:,second),[],4);
    halve(2).varValue(nElement,:) = nanvar(matrix(:,:,:,second),[],4);
    
    % Display progress
    if options.progress
        display_progress(iGroup, tGroup)
    end
    
end

% Elements of halve structure
halve(1).number = ceil(size(matrix,4)/2);
halve(2).number = floor(size(matrix,4)/2);

%% Set outputs

% Set output
varargout{1} = meanValue;
varargout{2} = stdValue;
varargout{3} = varValue;
varargout{4} = options;
varargout{5} = halve;

% Save data
if options.save
    save([options.outFolder '/' options.outName],...
        'meanValue','stdValue','varValue','options','halve');
end

% Display progress
if options.progress
    fprintf('Done\n');
end

end

