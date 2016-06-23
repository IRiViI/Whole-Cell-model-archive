% [track_clusters,cluster_list,archive] = find_clusters_archive(z,...
%   nClusters,sClusters,archive,strains,reactions_used,snapShotNumber)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 11/03/2016

function varargout = find_clusters_archive(varargin)

% Mandatory input arguments
z = varargin{1};                % dendogram information thing (z = linkage(y));
options.clusters = varargin{2}; % Number of clusters that should be found
options.minimal = varargin{3};  % The minimal number of objects in every cluster
archive = varargin{4};
options.included = varargin{5}; % The reactions/molecules included during the construction of the z matrix
options.snapShot = varargin{6}; % Number of the snap shop used in order to create the z matrix

% Default settings
options.sets = 1:length(archive.set);

% Adjust options
inputOptions = struct(varargin{7:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

% Find function
[trackClusters, clusterList, breakState] = find_clusters_function(z,options.minimal,options.clusters);

% Set output values if a solution is found
if breakState
    varargout{1} = trackClusters;
    varargout{2} = clusterList;
else
    fprintf('No clusters found with the current settings\n')
    varargout{1} = 0; varargout{2} = 0;
end

clusterOrder = 1;

% Add settings
archive.settings.clustering(clusterOrder).clusters = options.clusters;
archive.settings.clustering(clusterOrder).minimal = options.minimal;
archive.settings.clustering(clusterOrder).snapShot = options.snapShot;
archive.settings.clustering(clusterOrder).timePoint = ...
    archive.settings.snapShot(options.snapShot).timePoint;
archive.settings.clustering(clusterOrder).field = ...
    archive.settings.snapShot(options.snapShot).field;
archive.settings.clustering(clusterOrder).fileTypeTag = ...
    archive.settings.snapShot(options.snapShot).fileTypeTag;
archive.settings.clustering(clusterOrder).included = options.included;

% Update archive
cStrain = 1;
for iStrain = options.sets
    archive.set(iStrain).cluster(clusterOrder).value = [trackClusters(cStrain,:)];
    cStrain = cStrain + 1;
end

% Export archive
varargout{3} = archive;

end