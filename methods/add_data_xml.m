%Update the edges of a XML file for cytoscape. This scripts allows you to
%   set the thickness and the colours of the edges.
%   add_data_xml(xmlFile,reactionID,thicknessData,colorData).
%
%   xmlFile: File name of the xml file.
%   reactionID: List that contains all the names of the reactions. All the
%   eges of a specific reaction will get the same colour and thickness.
%   thicknessData: Data that's used to set the thickness of every egde of a
%   specific reaction. Values should be between 0 and 1.
%   colorData: Data for colouring the edges of a specific reaction. Values
%   should be between 0 and 1. If a value is "Nan" The colour will be set
%   to blue. A value of 0 will be set to red and values higher up to a
%   value of 1 become more green coloured.
%
%   'outputFile'    - Set the output file name. Default: Data_"inputFile".
%   'progress'      - Show progress. Default: 'On'
%   'coloring'      - Set the coloring of the lines. So far there are two
%                     options 'Default' and 'Jet'. (look at the code of
%                     this function to figures out how it works)
%   'metaboliteData'- include coloring information for the metabolites
%                     similar to thicknessData and colorData. Requires:
%                     'metaboliteID'.
%   'metaboliteID'  - The reactionID list

%   'clipped_metabolic_pathway_structure.xgmml'

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 02/12/2016

function add_data_xml(varargin)

% ----------- Process input arguments ----------- %

% Get options
if nargin > 4
    options = struct(varargin{5:end});
else
    options = struct();
end

% Get mandatory input arguments
options.inputFile = varargin{1};
options.reactID = varargin{2};
options.thicknessData = varargin{3};
options.colorData = varargin{4};

% Add default settings if no values are already given.
if ~isfield(options,'outputFile')
    options.outputFile = ['data_' options.inputFile];
end
if ~isfield(options,'coloring')
    options.coloring = 'Jet';
end
if ~isfield(options,'progress')
    options.progress = 'On';
end
if ~isfield(options,'metaboliteData')
    options.metabolite = 'Off';
else
    if isempty(options.metaboliteData)
        options.metabolite = 'Off';
    else
        options.metabolite = 'On';
    end
end
if ~isfield(options,'metabolicReaction')
    options.metabolicReaction = 'On';
end

options.compartments={'[c]','[]','[e]','[m]','[]','[]'};
% tc terminal organ cyt
% tm terminal organ mem
% d chromsome compt
% Test input arguments
if ~(length(options.colorData) == length(options.thicknessData))
    error('The list of color data and thickness data should be the same length');
end

% ----------- Process input arguments ----------- %

% Progress spacer
if strcmp(options.progress,'On')
    fprintf('Progress:\n')
end

% Load xml file
docNode = xmlread(options.inputFile);

% Get graph
graph = item(docNode,0);

% ----------- Reaction edges ----------- %

if strcmp(options.metabolicReaction,'On')
    
    % Find the source and target of every edge
    edgeInfo = getTypeOfTarget(graph,'edge',{'raw_source','raw_target'});
    
    % Create a single list in which both the points of an edge are saved
    edgePointList = {edgeInfo.raw_source; edgeInfo.raw_target}';
    
    % Total number of reactions
    tReaction = length(options.colorData);
    
    % Progress spacer
    if strcmp(options.progress,'On')
        fprintf('Add edge data\n     ')
    end
    
    % For every reaction
    for iReaction = 1:tReaction
        
        % Show progress
        if strcmp(options.progress,'On')
            display_progress(iReaction,tReaction)
        end
        
        % Color edge
        color = colorMapping(options.colorData(iReaction),options.coloring);
        
        % Thickness edge
        thickness = num2str(thicknessMapping(options.thicknessData(iReaction)));
        
        % Select the edges that belongs to the specific reaction
        edgeInList = find(sum(strcmp(options.reactID(iReaction),edgePointList),2));
        
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

% ----------- Metabolite nodes ----------- %
% (can be upgraded to handeling reactions aswell)

if strcmp(options.metabolite,'On')
    
    % Find the source and target of every edge
    nodeInfo = getTypeOfTarget(graph,'node','raw_label');
    
    % Total number of nodes
    t_node = length(nodeInfo);
    
    % Progress spacer
    if strcmp(options.progress,'On')
        fprintf('Add metabolite data:\n     ')
    end
    
    % For every node
    for i_node = 1:t_node
        
        % Show progress
        if strcmp(options.progress,'On')
            display_progress(i_node,t_node)
        end
        
        % Get node
        node = item(graph,nodeInfo(i_node).number);
        
        % Get graphics child
        graphicsInfo = getTypeOfTarget(node,'graphics','fill');
        graphics = item(node,graphicsInfo.number);
        
        % Look for the metabolite (if it is a metabolite)
        % Get the molecule number
        n_metabolite = find(strcmp(nodeInfo(i_node).raw_label(1:end-3),options.metaboliteID));
        if ~isempty(n_metabolite)
            
            % Get the compartment of the metabolite
            n_compartment = find(strcmp(nodeInfo(i_node).raw_label(end-2:end),options.compartments));
            
            % There is a change that a reaction slips trough the first
            % check, so check if it is indeed a compartement. Reactions do
            % not have a compartment
            if ~isempty(n_compartment)
                
                % Get color
                color = colorMapping(options.metaboliteData(n_metabolite,n_compartment),options.coloring);
                
                % Set color
                graphics.setAttribute('fill',color);
            else
                % It's probably a reaction
            end
        end
    end
    
end

% Save file
xmlwrite(options.outputFile,docNode);

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