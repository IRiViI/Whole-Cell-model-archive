function [trackClusters, clusterList, breakState] = find_clusters_function(z,minimal,clusters)
% z = z, minimal is minimal size of clusters, clusters is the number of
% clusters.

% initiate some variables
trackClusters = [];
breakState = false;

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
    numberOfClusters = sum(sortedList > minimal);
    
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
    if numberOfClusters == clusters
        breakState = true;
        break
    end
    
end

% Total number clusters
tTrackCluster = max(trackClusters(:,end));

% Number of splits
tSplit = size(trackClusters,2);

% Create cluster list
clusterList = [];
for iSplit = 1:tSplit
for iTrackCluster = 1:tTrackCluster
    clusterList(iTrackCluster,iSplit) = sum(trackClusters(:,iSplit)==iTrackCluster);
end
end

end