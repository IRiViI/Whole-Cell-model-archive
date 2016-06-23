% Add data values to WCM pathway visualization in Cytoscape.
%   3 input have to be set in order to change the value of a node or edge.
%   1) 'proteinComplexState': Set the state of the node or edge to true
%   2) 'proteinComplex': Give the ids of the elements of the nodes and 
%       edges
%   3) 'proteinComplexData': Give the values of the nodes and edges
%   On default: the 2nd input is set to the same id as the the state file 
%   outputs. So this options can be left untouched in some cases.
%   Expection: metabolicReaction has two different 'metabolicReactionData'
%   inputs instead 'colorData' and 'thicknessData'
%
%   Options:
%   Data input options:
%   'proteinComplexState'
%   'proteinComplex'
%   'proteinComplexData'
%   'proteinMonomerState'
%   'proteinMonomer'
%   'proteinMonomerData'
%   'metaboliteState'
%   'metabolite'
%   'metaboliteData'
%   'metabolicReactionState'
%   'metabolicReaction'
%   'colorData'
%   'thicknessData'
%   Other options:
%   'inputFile'                 - path to input xgmml file
%   'outputFile'                - name output file
%   'progress'                  - turn progress on or off (true or false)
%   
%   Example 1 (change color of protein nodes, of specific proteins):
%       archive.cyto_viz(...
%           'proteinComplexState',true,...
%           'proteinComplexData',proteinComplexData,...
%           'proteinComplex',proteinComplex);
%       % where:
%       proteinComplex(1).ID = 'DNA_GYRASE';    % Protein id
%       proteinComplex(1).form = 'mature';      % Protein form
%       proteinComplexData(1,:) = [1,1,1,1,1,1]; % Values for all 6
%                                                   locations
%       proteinComplex(2).ID = 'MG_215_TETRAMER';
%       proteinComplex(2).form = 'mature';
%       proteinComplexData(2,:) = [0,0,0,0,0,0];
%
%   Example 2 (Change color of all the metabolite nodes)
%       archive.cyto_viz(...
%           'metaboliteState',true,...
%           'metaboliteData',1-metaboliteData);
%       % 'metabolite' not required because it set on default to include
%       %   all the metabolites according to the state files
%       % where:
%       [stAveWTV,siAveWTV] = archive.average_value(1,2500,...
%           'Metabolite','counts',1:6);
%       metA = squeeze(stAveWTV);
%       metW = squeeze(siAveWTV);
%       metaboliteData = (metW(:,:,1) - metA(:,:))./metA(:,:);
%
%   Example 3 (Add metabolite and metabolic reaction data)
%       archive.cyto_viz('metaboliteState',true,...
%           'metaboliteData',1-squeeze(tmet),...
%           'metabolicReactionState',true,...
%           'colorData',1-squeeze(tfluxs))
%       where:
%       [stAveWTV,siAveWTV] = archive.average_value(1,2500,...
%           'Metabolite','counts',1:6);
%       [stAveWT,siAveWT] = archive.average_value(1,2500,...
%           'MetabolicReaction','fluxs',1);
%       metA = squeeze(stAveWTV);
%       metW = squeeze(siAveWTV);
%       tmet = (metW(:,:,1) - metA(:,:))./metA(:,:);
%       tfluxs = (siAveWT(:,:,:,1) - stAveWT)./stAveWT;

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/12/2716

function add_data_xml(archive, varargin)
%% ----------- Process input arguments ----------- %
%% Default settings
options.inputFile = 'Default_layout.xgmml';
options.outputDir = [archive.settings.dir '/' 'output'];
for iReaction = 1:length(archive.info.MetabolicReaction)
options.metabolicReaction(iReaction).ID = archive.info.MetabolicReaction(iReaction).ID;
end
for iMetabolite = 1:length(archive.info.Metabolite)
options.metabolite(iMetabolite).ID = archive.info.Metabolite(iMetabolite).ID;
end
for iProtein = 1:length(archive.info.ProteinComplex)
options.proteinComplex(iProtein).ID = archive.info.ProteinComplex(iProtein).ID;
options.proteinComplex(iProtein).form = archive.info.ProteinComplex(iProtein).form;
end
for iProtein = 1:length(archive.info.ProteinMonomer)
options.proteinMonomer(iProtein).ID = archive.info.ProteinMonomer(iProtein).ID;
options.proteinMonomer(iProtein).form = archive.info.ProteinMonomer(iProtein).form;
end
options.proteinMonomerID = {archive.info.ProteinMonomer.ID};
options.node = true;
options.progress = true;
options.metaboliteState = false;
options.metabolicReactionState = false;
options.proteinComplexState = false;
options.proteinMonomerState = false;
options.coloring = 'Jet';
options.compartments={'[c]','[]','[e]','[m]','[]','[]'};
% tc terminal organ cyt
% tm terminal organ mem
% d chromsome compt
tReaction = length(archive.info(1).MetabolicReaction);
options.colorData = 0.5*ones(1,tReaction);
options.thicknessData = 0.5*ones(1,tReaction);

% Apply flux values
if options.metabolicReactionState
    if isfield(options,'colorData')
        options.colorData = options.colorData;
    end
    if isfield(options,'thicknessData')
        options.thicknessData = options.thicknessData;
    end
end

%% Adjust options

inputOptions = struct(varargin{:});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

%% ----------- Process settings ----------- %

% Test input arguments
if ~(length(options.colorData) == length(options.thicknessData))
    error('The list of color data and thickness data should be the same length');
end

% Create output file name
if ~isfield(options,'outputFile')
tmp = strsplit_archive(options.inputFile,'/');
options.outputFile = [options.outputDir '/' 'data_' tmp{end}];
end

% Add library
addpath([archive.settings.dir '/' 'functions/others/hex_and_rgb_v2']);

%% ----------- Process input arguments ----------- %

% Progress spacer
if options.progress
    fprintf('Progress:\n')
end

% Load xml file
docNode = xmlread(options.inputFile);

% Get graph
graph = item(docNode,0);

%% ----------- Edges ----------- %

graph = add_reaction_edge_info(graph,options);

%% ----------- Nodes ----------- %
% (can be upgraded to handeling reactions aswell)

graph = add_node_info(graph,options);

%% ----------- Save file ----------- %

xmlwrite(options.outputFile,docNode);

% Progress
fprintf(['Done: ' options.outputFile '\n'])

end

function color = colorMapping(value,colorScheme)
% Determine the right color. Value should be between 0 and 1

switch colorScheme
    case 'Default'
        try
            color = rgb2hex([1-value value 0]);
        catch
            color = rgb2hex([0 0 1]);
        end
    case 'Jet'
        matrix = jet;
        try
            % Set 0 to 1 and 1 to 64. Where 1 and 64 is the range of jet.
            color = rgb2hex(matrix(floor(value*(length(matrix)-1))+1,:));
        catch % Exceptions
            
            if value < 0 && ~(value == -Inf)
                color = rgb2hex([0 0 0.5]);
            elseif value > 1 && ~(value == Inf)
                color = rgb2hex([0.5 0 0]);
            elseif value == Inf
                color = rgb2hex([0 0 0.5]);
            elseif value == -Inf
                color = rgb2hex([0.5 0 0]);
            elseif isnan(value)
                color = rgb2hex([1 1 1]);
            else
                %                 color = rgb2hex([0 0 0]);
                error('non valid value')
            end
        end
end

end

function thickness = thicknessMapping(value)
% Determine the right thickness. Value should be between 0 and 1
if value >= 0 && value <= 1
    thickness = num2str(value*9+1);
elseif value < 0
    thickness = 1;
elseif value > 1
    thickness = 10;
elseif value == inf
    thickness = 10;
elseif value == -inf;
    1;
else
    thickness = 0;
end
end

function graph = add_reaction_edge_info(graph,options)

if options.metabolicReactionState
    
    % Find the source and target of every edge
    edgeInfo = getTypeOfTarget(graph,'edge',{'raw_source','raw_target'});
    
    % Create a single list in which both the points of an edge are saved
    edgePointList = {edgeInfo.raw_source; edgeInfo.raw_target}';
    
    % Total number of reactions
    tReaction = length(options.colorData);
    
    % Progress spacer
    if options.progress
        fprintf('Add edge data\n     ')
    end
    
    % For every reaction
    for iReaction = 1:tReaction
        
        % Show progress
        if options.progress
            display_progress(iReaction,tReaction)
        end
        
        % Color edge
        color = colorMapping(options.colorData(iReaction),options.coloring);
        
        % Thickness edge
        thickness = num2str(thicknessMapping(options.thicknessData(iReaction)));
        
        % Select the edges that belongs to the specific reaction
        edgeInList = find(sum(strcmp(options.metabolicReaction(iReaction).ID,edgePointList),2));
        
        % Get the numbers of the edges belonging to the xml file
        edgeNumbers = [edgeInfo(edgeInList).number];
        
        % Number of edges
        tEdge = length(edgeNumbers);
        
        % For every edge
        for iEdge = 1:tEdge
            
            % Get edge
            edge = item(graph,edgeNumbers(iEdge));
            
            % Get graphics child
            graphicsInfo = getTypeOfTarget(edge,'graphics','fill');
            graphics = item(edge,graphicsInfo.number);
            
            % Set color
            graphics.setAttribute('fill',color);
            
            % Set thinckness
            graphics.setAttribute('width',thickness);
            
        end
        
    end
    
end

end

function graph = add_node_info(graph,options)

if options.node
    
    % Find the source and target of every edge
    nodeInfo = getTypeOfTarget(graph,'node','raw_id');
    
    % Total number of nodes
    t_node = length(nodeInfo);
    
    % Progress spacer
    if options.progress
        fprintf('Add node data:\n     ')
    end
    
    % For every node
    for i_node = 1:t_node
        
        % Show progress
        if options.progress
            display_progress(i_node,t_node)
        end
        
        % Get node
        node = item(graph,nodeInfo(i_node).number);
        
        % Get graphics child
        graphicsInfo = getTypeOfTarget(node,'graphics','fill');
        graphics = item(node,graphicsInfo.number);
        
        % Look for the molecule
        n_molecule = [];
        n_data = '';
        n_compartment = [];
        if options.metaboliteState
        % Get the molecule number
        n_molecule = find(strcmp(nodeInfo(i_node).raw_id(1:end-3),{options.metabolite.ID}));
        n_data = 'metaboliteData';
        end
        % Check if it's a protein Complex if it's not a metabolite
        if isempty(n_molecule) && options.proteinComplexState
            n_molecule = find(strcmp(nodeInfo(i_node).raw_id(1:end-3),{options.proteinComplex.ID}));
            % Get the mature protein
            n_form = strcmpi({options.proteinComplex(n_molecule).form},'mature');
            n_molecule = n_molecule(n_form);
            if ~isempty(n_molecule)
            n_data = 'proteinComplexData';
            end
        end
        % Check if it's a protein Monomer if it's not a metabolite or
        % protein Complex
        if isempty(n_molecule) && options.proteinMonomerState
            n_molecule = find(strcmp(nodeInfo(i_node).raw_id(1:end-3),{options.proteinMonomer.ID}));
            % Get the mature protein
            n_form = strcmpi({options.proteinMonomer(n_molecule).form},'mature');
            n_molecule = n_molecule(n_form);
            if ~isempty(n_molecule)
            n_data = 'proteinMonomerData';
            end
        end
            
        if ~isempty(n_molecule)
            
            % Get the compartment of the molecule
            n_compartment = find(strcmp(nodeInfo(i_node).raw_id(end-2:end),options.compartments));
            
            % There is a change that a reaction slips trough the first
            % check, so check if it is indeed a compartement. Reactions do
            % not have a compartment.
            if ~isempty(n_compartment)
                
                % Get color
                color = colorMapping(options.(n_data)(n_molecule,n_compartment),options.coloring);
                
                % Set color
                graphics.setAttribute('fill',color);
            else
                % It's probably a reaction
            end
        end
    end
    
end

end