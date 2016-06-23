% Signal to noise is determined by comparing the values between different
% simulations. It might be better to use the standard deviation of single
% simulations instead.
%
%   Example:
%   clustering_extensive_archive(archive,[1,2,3,4],1,'newFigure',false,...
%       'clustering',true,'clusterNumber',[4,3,2,3],...
%       'clusterSize',[15,15,30,15],'combineCluster',true,...
%       'combinedClusterNumber',5,'combinedClusterSize',10);

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/27/2016

function varargout = clustering_extensive_archive(archive, snapShot, reference, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.snapShot = snapShot;
options.reference = reference;
options.SNR = 2;            % The SNR input can be both a scalar or a vector with the same length as the number of snapshots
options.compartments = [true,true,false,true,true,true]; % Do not include extracellular molecules
options.tresshold = 'Default';
options.dendrogram = true;          % Start with analyzing the hiearchical clusters
options.clustering = false;         % Find the clusters
options.minClusterSize = 5;
options.clusterNumber = [];
options.clusterSize = [];
options.newFigure = true;
options.combineCluster = false;

% Adjust settings
tmpOptions = struct(varargin{1:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

%% Process settings

% Total number of snapshots to analyze
options.tSnapShot = length(options.snapShot);

% Get all sets in case that all sets should be analyzed
if options.set == 0
    options.set = 1:archive.sets;
end

% Process signal to noise settings
options = process_SNR_settings(options);

% Tress hold values
options = process_tresshold_settings(archive,options);

% Save options
out.options = options;

%% Run code

% Display dendrograms for identifying clusters
if options.dendrogram
    [DG,out] = create_dendogram(archive,out,options);
end

% Process sub clustering
if options.clustering
    out = process_sub_clustering(archive, out,DG,options);
end

% Combine clusters
if options.combineCluster
    % Display progress
    if options.progress
        fprintf('Combine clusters\n');
    end
    [archive, out] = combine_clusters(archive, out, options);
end

% Display progress
if options.progress
    fprintf('Done\n');
end

%% Set outputs

% Set output
varargout{1} = out;
varargout{2} = archive;

end

function options = process_SNR_settings(options)
% Every snap should should have its own signal-to-noise ratio value. This
% function changes a scalar value to a vector with the same number of
% elements as the number of snap shots.

if isscalar(options.SNR)
    options.SNR = options.SNR .* ones(1,length(options.snapShot));
elseif length(options.SNR) == length(options.snapShot)
    % Do nothing, the SNR values already match.
else
    error('Number of signal-to-noise ratio values does not match the number of snap shots');
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
            snapShot = options.snapShot(iSnapShot);
            field = archive.settings.snapShot(snapShot).field;
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
    options.tresshold = options.tresshold .* ones(1,length(options.snapShot));
elseif length(tresshold) == length(options.snapShot) && ~ischar(tresshold);
    options.tresshold = tresshold;
else
    error('Number of tresshold values does not match the number of snap shots');
end
end

function matrix = apply_filters(archive,iSnapShot,options)
% Create a matrix with the same size as the snap shot matrix in which
% every element of the matrix indicates true if the molecule or reaction
% meets the filter criteria

% ------------ Preprocessing

% Current snapShot number
snapShot = options.snapShot(iSnapShot);

% Get signal and noise values of snap shot
signalRef = archive.set(options.reference).snapShot(snapShot).averageValues;
noiseRef = archive.set(options.reference).snapShot(snapShot).stdValues;

% Dimensions matrix
[x,~] = size(signalRef);

% ------------ Create filters

% Apply SNR criteria
filter(:,:,1) = (signalRef./noiseRef) > options.SNR(iSnapShot);

% Apply compartment criteria
if size(filter,2) == 6;
    filter(1:x,options.compartments,end+1) = true;
end

% ------------ Post filters processing

% Combine filters
matrix = prod(filter,3);

end

function clusterSizes = get_cluster_sizes(clusters)
% Get the cluster sizes at all cluster split junkions
tJunktion = size(clusters,2); % Number of split junkions in clustering
for iJunktion = 1:tJunktion
    tCluster = max(clusters(:,iJunktion));
    for iCluster = 1:tCluster
        clusterSizes(iCluster,iJunktion) = sum(clusters(:,iJunktion)==iCluster);
    end
end
end

function out = process_sub_clustering(archive, out,DG, options)

for iSnapShot = 1:options.tSnapShot
    
    % Get snapshot number
    snapShot = options.snapShot(iSnapShot);
    
    % Make output structure
    out.snapShot(iSnapShot).tresshold = options.tresshold(iSnapShot);
    out.snapShot(iSnapShot).SNR = options.SNR(iSnapShot);
    
    % Progress
    if options.progress
        fprintf(['Clustering snap shot ' num2str(snapShot) '\n']);
    end
    
    % Find clusters
    [~,clusters,~] = archive.clustering(DG(iSnapShot).z,...
        options.clusterNumber(iSnapShot),...
        options.clusterSize(iSnapShot),...
        options.set,...
        snapShot);
    
    % Get cluster sizes at all the split junktions
    clusterSizes = get_cluster_sizes(clusters);
    
    % Save outputs
    out.snapShot(iSnapShot).clusterSizes = clusterSizes;
    out.snapShot(iSnapShot).clusterHistory = clusters;
    
    % Number of junktions
    tJunktion = size(clusterSizes,2);
    
    for iJunktion = 1:tJunktion
        
        
        % Put sets in clusters
        % List of clusters that meet criteria
        lCluster = find(clusterSizes(:,iJunktion) > options.minClusterSize);
        % Total number of clusters
        tCluster = size(lCluster,1);
        
        % Initate fields for all sets
        for iSet = 1:archive.sets
            out.snapShot(iSnapShot).set(iSet).cluster(1:tCluster,iJunktion) = false;
        end
        
        % assign sets to cluters
        for iCluster = 1:tCluster;
            
            % Find sets that agree with the cluster number
            sets = options.set(clusters(:,iJunktion) == lCluster(iCluster));
            
            % Set values
            for iSet = sets
                out.snapShot(iSnapShot).set(iSet).cluster(iCluster,iJunktion) = true;
            end
        end
        
    end
    
end
end

function [DG,out] = create_dendogram(archive,out,options)
for iSnapShot = 1:options.tSnapShot
    
    % Get snapshot number
    snapShot = options.snapShot(iSnapShot);
    
    % Filter data
    filter = apply_filters(archive,iSnapShot,options);
    
    % Set values to output
    out.snapShot(iSnapShot).number = snapShot;
    out.snapShot(iSnapShot).field = archive.settings.snapShot(snapShot).field;
    out.snapShot(iSnapShot).timePoint = archive.settings.snapShot(snapShot).timePoint;
    out.snapShot(iSnapShot).fileTypeTag = archive.settings.snapShot(snapShot).fileTypeTag;
    
    % New figure
    if options.newFigure
        figure;
    end
    
    % Create dendogram
    DG(iSnapShot).z = archive.dendrogram(snapShot,...
        'include',filter,...
        'compare','average',...
        'selection',options.tresshold(iSnapShot),...
        'set', options.set);
    
    % figure options.
    set(gcf,'Position',[400,300,400,200])
    title(['Dendogram of snap shot: ' num2str(snapShot) ': ' out.snapShot(iSnapShot).field]);
    
end
end

function [archive, out] = combine_clusters(archive, out, options)
% Add: combinedClusterNumber, combinedClusterSize

% Combine the sub clusters
for iSet = options.set
    for iSnapShot = 1:options.tSnapShot
        values = out.snapShot(iSnapShot).set(iSet).cluster;
        matrix(iSnapShot).values = values(:)./size(values,2);
    end
    out.set(iSet).cluster = cat(1,matrix.values);
end
X = double(cat(2,out.set.cluster))';

% Create z structure
Y = pdist(X);
Z = linkage(Y);

% Create dendogram
if options.newFigure
    figure
end
D = dendrogram(Z,0);

% Find clusters
[trackClusters,~, breakState] = find_clusters_function(Z,...
    options.combinedClusterSize,...
    options.combinedClusterNumber);
if breakState
else
    fprintf('No clusters found with the current settings\n');
end

% Total number clusters
tTrackCluster = max(trackClusters(:,end));

% Create cluster list
for iTrackCluster = 1:tTrackCluster
    clusterList(iTrackCluster) = sum(trackClusters(:,end)==iTrackCluster);
end

out.trackClusters = trackClusters;
out.clusterList = clusterList;

lCluster = find((clusterList > options.minClusterSize));
tCluster = length(lCluster);

% Initiate combined cluster values
for iSet = 1:archive.sets
    out.set(iSet).combinedCluster = 0;
     archive.set(iSet).combinedCluster = 0;
end

for iCluster = 1:tCluster
    cluster = lCluster(iCluster);
    % list of sets with the cluster number
    lSet = find(trackClusters(:,end)==cluster)';
    for iSet = lSet
        out.set(iSet).combinedCluster = iCluster;
        archive.set(iSet).combinedCluster = iCluster;
    end
end

end