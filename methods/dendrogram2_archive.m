% Improved version of dendrogram_archive.
%   [z,D,options] = dendrogram2_archive(archive, snapShot, refSet,...)
%   snapShot can be an array of snap shot numbers
%   referense is the referense set to use filtering
%   z is z matrix of dendrogram
%   D is the output of the dendrogram
%   options contains the settings.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/31/2016

function varargout = dendrogram2_archive(archive, snapShot, refSet, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.target = 'set';     % Default: cluster sets
options.lSnapShot = snapShot;
options.tSnapShot = length(options.lSnapShot); % Don't change!
options.referenceSet = refSet;
options = set_field_filters(archive,options); % Snap shot specific filters
options.tresshold = 1;
options.SNR = 2;
options.factor = 1;
options.dargin = {};            % arguments for dendrogram

% Change settings
tmpOptions = struct(varargin{1:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

%% Process settings

% Determine which sets of the archive should be analyzed
if options.set == 0
    options.tSet = archive.sets;    % Number of sets to analyze
    options.lSet = 1:options.tSet;          % List of sets to analyze
else
    options.lSet = options.set;     % List of sets to analyze
    options.tSet = length(options.lSet);    % Number of sets to analyze
end

% Determine number of simulations to analyze (required for progess report)
if options.simulation == 0
    options.ttSimulation = archive.simulations(options.lSet);
else
    options.ttSimulation = tSet * length(options.simulation);
end

% Set filter values
options = process_tresshold_settings(archive,options);
options = process_target_settings(options,'SNR');

% Set factor values
options = process_target_settings(options,'factor');

%% Main process

% -------- Filters -------- %

% Create filters
options = create_filters(archive,options);

% -------- Data -------- %

% Create reference vector (used for filtering
[refVector,snapShotVector] = create_ref_vector(archive, options);
options.refVector = refVector;
options.snapShotVector = snapShotVector;

% Create matrix which should be used for clustering
dataMatrix = create_data_matrix(archive,options);
options.dataMatrix = dataMatrix;

% Reference matrix
dSize = size(dataMatrix);
refMatrix = refVector*ones(1,dSize(2));

% -------- Scale values -------- %

% Compare values with reference
data = (dataMatrix-refMatrix)./refMatrix;

% Rescale values
x = log(abs(data)+1).*((data>0)*2-1);

% Adjust for snap shot type
for iSnapShot = 1:options.tSnapShot
    x(snapShotVector==iSnapShot,:) = x(snapShotVector==iSnapShot,:) * options.factor(iSnapShot);
end

if strcmp(options.progress,'On')
    fprintf('Done\n')
end

% -------- Create dendrogram -------- %
    
% Prepare dendrogram
x = x';
y = pdist(x);
z = linkage(y, 'complete');
D = dendrogram(z,0,options.dargin{:});

% Aesthetics
set(gcf,'Position',[400,300,600,200])
title([archive.settings.snapShot(options.lSnapShot).field ...
    ', at ' num2str(archive.settings.snapShot(options.lSnapShot).timePoint) ' s'])
set(gca,'XTickLabel','')

% Display progress
if options.progress
    fprintf('Done\n');
end

%% Set outputs

% Set output
varargout{1} = z;
varargout{2} = D;
varargout{3} = options;

end

function options = set_field_filters(archive, options)

% initiate filter structure
options.filter.snapShot = [];

% Make a filter for every snap shot
for iSnapShot = 1:options.tSnapShot
    
    
    nSnapShot = options.lSnapShot(iSnapShot);
    
    % Field of snap shot
    field = archive.settings.snapShot(nSnapShot).field;
    
    % Size of snap shot
    mSize = size(archive.set(options.referenceSet).snapShot(nSnapShot).averageValues);
    
    % Initiate filter matrix
    matrix = ones(mSize);
    
    % Adjust filter
    if strcmpi(field,'MetabolicReaction')
        % Include everything
    elseif strcmpi(field,'Metabolite')
        % Exclude extracellular
        matrix(:,3) = 0;
    elseif strcmpi(field,'ProteinComplex')
        % Exclude extracellular
        matrix(:,3) = 0;
    elseif strcmpi(field,'ProteinMonomer')
        % Exclude extracellular
        matrix(:,3) = 0;
    elseif strcmpi(field,'Rna')
        matrix(:,3) = 0;
    else
        error(['snap shot field ' num2str(nSnapShot) ' does not have a filter']);
    end
    
    % Save filter
    options.filter.snapShot(iSnapShot).field = boolean(matrix);
    
end

end

function matrix = create_data_matrix(archive, options)
% Create one single matrix with all the data values that should be used for
% the clustering

% Analyze sets
for iSet = options.lSet
    
    % Clear new Data structure
    newData = [];
    
    % Go through all targeted snap shots
    for iSnapShot = 1:options.tSnapShot
        
        nSnapShot = options.lSnapShot(iSnapShot);
        
        % field of snap shot
        field = archive.settings.snapShot(nSnapShot).field;
        
        % Get filter
        filter = options.filter.snapShot(iSnapShot).total;
      
        % Add snap shot data to the set structure
        newData(end+1).values = archive.set(iSet).snapShot(nSnapShot).averageValues(filter);
        
    end
    
    % Combine snapShots
    set(iSet).data = cat(1,newData.values);
    
end

% Create one big matrix
matrix = cat(2,set.data);

end

function [refVector,snapShotVector] = create_ref_vector(archive,options)
% Create a vecot of the reference set

for iSnapShot = 1:options.tSnapShot
   
    % Current snapshot
    nSnapShot = options.lSnapShot(iSnapShot);
    
    % Get filter
    filter = options.filter.snapShot(iSnapShot).total;
    
    % Get values
    vector(iSnapShot).values = archive.set(options.referenceSet).snapShot(nSnapShot).averageValues(filter);
    
    % Track the snap shot of the elements
    vector(iSnapShot).field = iSnapShot * ones(size(vector(iSnapShot).values));
    
end

% Combine the values
refVector = cat(1,vector.values);
snapShotVector = cat(1,vector.field);

end

function options = create_filters(archive,options)
% Create the filter that should be applied on the data
% Use the reference strain and field specific filters.

% Process filter for all snap shots
for iSnapShot = 1:options.tSnapShot
    
    % Current snap shot
    nSnapShot = options.lSnapShot(iSnapShot);
    
    % Get values of snap shot 
    values = archive.set(options.referenceSet).snapShot(nSnapShot).averageValues;
    std = archive.set(options.referenceSet).snapShot(nSnapShot).stdValues;
    
    % Size of snap shot
    mSize = size(values);
    
    % Initiate filter matrix
    matrix = ones(mSize);
    
    % Tresshold filter
    tmpFilter = values > options.tresshold(iSnapShot);
    matrix(~tmpFilter) = 0;
    options.filter.snapShot(iSnapShot).tresshold = boolean(tmpFilter);
    
    % SNR filter
    tmpFilter = values./std > options.SNR(iSnapShot);
    matrix(~tmpFilter) = 0;
    options.filter.snapShot(iSnapShot).SNR = boolean(tmpFilter);
    
    % field specific filter
    tmpFilter = options.filter.snapShot(iSnapShot).field;
    matrix(~tmpFilter) = 0;
    
    % Final filter
    options.filter.snapShot(iSnapShot).total = boolean(matrix);
    
end


end

function options = process_target_settings(options,target)
% Every snap should should have its own signal-to-noise ratio value. This
% function changes a scalar value to a vector with the same number of
% elements as the number of snap shots.

if isscalar(options.(target))
    options.(target) = options.(target) .* ones(1,options.tSnapShot);
elseif length(options.(target)) == options.tSnapShot
    % Do nothing, the SNR values already match.
else
    error(['Number of ' target ' values does not match the number of snap shots']);
end
end

function options = process_tresshold_settings(archive,options)
% Set a tresshold value for every snap shot

tresshold = options.tresshold;
options.tresshold = [];

if ischar(tresshold)
    % Default settings
    if strcmp(tresshold,'Default')
        for iSnapShot = 1:options.tSnapShot
            nSnapShot = options.snapShot(iSnapShot);
            field = archive.settings.snapShot(nSnapShot).field;
            if strcmpi(field,'MetabolicReaction')
                options.tresshold(iSnapShot) = 1;
            elseif strcmpi(field,'Metabolite')
                options.tresshold(iSnapShot) = 1;
            elseif strcmpi(field,'ProteinComplex')
                options.tresshold(iSnapShot) = 1;
            elseif strcmpi(field,'ProteinMonomer')
                options.tresshold(iSnapShot) = 1;
            elseif strcmpi(field,'Rna')
                options.tresshold(iSnapShot) = 1;
            else
                error(['snap shot field ' num2str(iSnapShot) ' not recognized']);
            end
        end
    end
    % Process manually set tresshold values
    % Match the number of vector elements to the number of snap shots
elseif isscalar(tresshold) && ~ischar(tresshold)
    options.tresshold = tresshold .* ones(1,options.tSnapShot);
elseif length(tresshold) == options.tSnapShot && ~ischar(tresshold);
    options.tresshold = tresshold;
else
    error('Number of tresshold values does not match the number of snap shots');
end
end