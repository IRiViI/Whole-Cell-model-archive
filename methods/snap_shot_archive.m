% Make a snap shot of the simulations at a certain time point
%
%   archive = snap_shot(archive,time_point,field,subfield)
%
%   'fileTypeTag'   - Make a snap shot of a certain file. Default:
%                     'mean_compressed'
% Output: (element, time, compartment)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/05/2016

function archive = snap_shot_archive(varargin)

% Get mandatory inputs
archive = varargin{1};
options.timePoint = varargin{2};
options.field = varargin{3};
options.subfield = varargin{4};

% Get optional inputs
if nargin > 4
    input_options = struct(varargin{5:end});
else
    input_options = struct();
end

% Set the other options
if isfield(input_options,'fileTypeTag')
    options.fileTypeTag = input_options.fileTypeTag;
else
    options.fileTypeTag = 'mean_compressed_state.mat';
end
if isfield(input_options,'compartment')
    options.compartment = input_options.compartment;
else
    options.compartment = 1;
end

% Get whether the fileTypeTag is direct or not
C = strsplit(options.fileTypeTag,'.');
if strcmp(C(end),'mat')
   options.direct = 'On'; 
else
    options.direct = 'Off'; 
end

% Number of sets
tStrain = length(archive.set);

% Progress spacer
fprintf('Progress:\n     ');

% Add settings to archive (make entry first if required)
if isfield(archive.settings,'snapShot')
    snapShotNumber = length(archive.settings.snapShot)+1;
else
    snapShotNumber = 1;
end

archive.settings.snapShot(snapShotNumber).timePoint = options.timePoint;
archive.settings.snapShot(snapShotNumber).field = options.field;
archive.settings.snapShot(snapShotNumber).fileTypeTag = options.fileTypeTag;

% For every set
for iStrain = 1:tStrain
    
    % Display progress
    display_progress(iStrain,tStrain)
    
    try
        % Get the values at specific time point of correct file
        [~,fullMatrix] = average_value_archive(archive,iStrain,options.timePoint,...
            options.field,options.subfield,options.compartment,...
            'fileTypeTag',options.fileTypeTag,...
             'direct',options.direct);
        
        % Number of simulations in set
        tSimulation = length(archive.set(iStrain).simulation);
        
        % For every simulation
        for iSimulation = 1:tSimulation
            % Set the values
            archive.set(iStrain).simulation(iSimulation).snapShot(snapShotNumber).values =...
                fullMatrix(:,:,:,iSimulation);
        end
        
        % Set average values of the set
        archive.set(iStrain).snapShot(snapShotNumber).averageValues = nanmean(fullMatrix,4);
        archive.set(iStrain).snapShot(snapShotNumber).stdValues = nanstd(fullMatrix,[],4);
    catch
        % Give warning when no values are found
        warningMessage = sprintf('It''s likely that set %d does not have the correct files',iStrain);
        warning(warningMessage);
    end
    
end

end