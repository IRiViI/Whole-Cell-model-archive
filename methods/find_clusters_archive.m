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

% Settings for adding the clustering information to the archive
if nargin > 3
    options.updatearchive = 1;
    if nargin > 3
        archive = varargin{4};
        options.strains = 1:length(archive.set);
        options.included = varargin{5}; % The reactions/molecules included during the construction of the z matrix
        options.snapShot = varargin{6}; % Number of the snap shop used in order to create the z matrix
    else
        error('Not enough input arguments for updating the achive')
    end
else
    options.updatearchive = 0;
end

% initiate some variables
breakState = 0;
trackClusters = [];

% Test different number of clusters
for iCluster = 1:length(z) % Start from 1 because of tracking cluster formation
    
    % Clustering
    t = cluster(z,'maxclust',iCluster);
    
    % Check the size of every cluster
    for iCheck = 1 : iCluster
        clusterList(iCheck) = sum(t == iCheck);
    end
    
    % Sort the cluster on their size (size matters here)
    sortedList = sort(clusterList);
    
    % The number enough large clusters
    numberOfClusters = sum(sortedList > options.minimal);
    
    % Check if there is an additonal cluster found
    try
        if numberOfClusters > previousNumberOfClusters
            % Add new cluster state to cluster track clusters list
            trackClusters(:,end+1) = t;
        end
    catch
    end
    % Save information for next time
    previousNumberOfClusters = numberOfClusters;
    
    % Check if there are enough large-enough-clusters found
    if numberOfClusters == options.clusters
        breakState = 1;
        break
    end
    
end

% Set output values if a solution is found
if breakState
    varargout{1} = trackClusters;
    varargout{2} = clusterList;
else
    fprintf('No clusters found with the current settings\n')
    varargout{1} = 0; varargout{2} = 0;
end

% Update the archive structure
if options.updatearchive
    
%     % Check now manyth clustering this is
%     if isfield(archive.settings,'clustering')
%         clusterOrder = length(archive.settings.clustering);
%     else
%         clusterOrder = 1;
%     end

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
    for iStrain = options.strains
        archive.set(iStrain).cluster(clusterOrder).value = [trackClusters(cStrain,:)];
        cStrain = cStrain + 1;
    end
    
    % Export archive
    varargout{3} = archive;
end

end