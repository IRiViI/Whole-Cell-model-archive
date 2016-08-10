function options = cytoscape_structure_archive(archive,varargin)

%% Settings

% Default settings
options.outName = 'output/cytoscape_WCM_pathway';
options.id = 'tXML';
options.name = 'test_XML_making';
options.metaid = 'geen idee';
options.layout = [];

% Default graphics settings molecule nodes
moleculeNode.h = 40;
moleculeNode.w = 40;
moleculeNode.x = 'shift_80';
moleculeNode.y = 80;
moleculeNode.type = 'ELLIPSE';
moleculeNode.fill = '#ffffff';
moleculeNode.width = 5;
moleculeNode.outline = '#ff0000';
options.moleculeNode = moleculeNode;

% Default graphics settings reaction nodes
reactionNode.h = 40;
reactionNode.w = 40;
reactionNode.x = 'shift_80';
reactionNode.y = 80;
reactionNode.type = 'DIAMOND';
reactionNode.fill = '#cccccc';
reactionNode.width = 5;
reactionNode.outline = '#c0c0c0';
options.reactionNode = reactionNode;

% Default graphics settings molecule nodes
proteinNode.h = 40;
proteinNode.w = 40;
proteinNode.x = 'shift_80';
proteinNode.y = 80;
proteinNode.type = 'ELLIPSE';
proteinNode.fill = '#ffffff';
proteinNode.width = 5;
proteinNode.outline = '#06A4A4';
options.proteinNode = proteinNode;

% Default graphics settings reaction edges
reactionEdges.width = 1;
reactionEdges.fill = '#000000';
options.reactionEdges = reactionEdges;

% Default graphics settings protein edges
proteinEdges.width = 1;
proteinEdges.fill = '#06A4A4';
options.proteinEdges = proteinEdges;

% Adjust options
inputOptions = struct(varargin{1:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

%% Preprocessing

% process layout structure (change format)
if ~isempty(options.layout)
    options.layout = format_layout(options.layout);
end

%% Initiation process

% Initiate xml structure
docNode = com.mathworks.xml.XMLUtils.createDocument('graph');
docNode = initial_xml_settings(docNode);

% Get graph structure
graph = item(docNode,0);

%% Create nodes

% Progress
fprintf('Add molecules nodes to xml file:\n     ');

% Add molecule nodes
graph = add_molecule_nodes(docNode, graph, archive, options);

% Progress
fprintf('Add reaction nodes to xml file:\n     ');

% Add reaction nodes
graph = add_reaction_nodes(docNode, graph, archive, options);

% Progress
fprintf('Add protein nodes to xml file:\n     ');

% Add protein nodes
graph = add_protein_nodes(docNode, graph, archive, options);

%% Process nodes

% Make a list of all the current nodes and find their location in the
% schematic. This will help to make the edges between the diffferent nodes.
lNode = process_nodes(graph);

%% Create Edges

% Progress
fprintf('Add reaction edges to xml file:\n     ');

% Add reaction edges
graph = add_reaction_edges(docNode, graph, archive, options, lNode);

% Progress
fprintf('Add protein-reaction edges to xml file:\n     ');

% Add reaction edges
graph = add_protein_edges(docNode, graph, archive, options, lNode);

% Save sml structure
xmlwrite([options.outName '.xgmml'],docNode);

% Progress
fprintf(['Done: ' options.outName '.xgmml' '\n'])

end

function docNode = initial_xml_settings(docNode)

% Create graph element (No idea what most this do but they are important!)
graph = docNode.getDocumentElement;
graph.setAttribute('label','Whole_Cell_Model');
graph.setAttribute('xmlns:xlink','http://www.w3.org/1999/xlink');
graph.setAttribute('xmlns:cy','http://www.cytoscape.org');
graph.setAttribute('xmlns:dc','http://purl.org/dc/elements/1.1/');
graph.setAttribute('xmlns:rdf','http://www.w3.org/1999/02/22-rdf-syntax-ns#');
graph.setAttribute('xmlns','http://www.cs.rpi.edu/XGMML');
graph.setAttribute('directed','1');

% Create att elements (because why not ay)
att = docNode.createElement('att');
att.setAttribute('name','documentVersion');
att.setAttribute('value','1.1');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','string');
att.setAttribute('name','backgroundColor');
att.setAttribute('value','#ffffff');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','real');
att.setAttribute('name','GRAPH_VIEW_ZOOM');
att.setAttribute('value','1');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','real');
att.setAttribute('name','GRAPH_VIEW_CENTER_X');
att.setAttribute('value','80');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','real');
att.setAttribute('name','GRAPH_VIEW_CENTER_Y');
att.setAttribute('value','0');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','boolean');
att.setAttribute('name','NODE_SIZE_LOCKED');
att.setAttribute('value','true');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','string');
att.setAttribute('name','sbml id');
att.setAttribute('value','tXML');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','string');
att.setAttribute('name','sbml metaId');
att.setAttribute('value','WC_Model');
graph.appendChild(att);

att = docNode.createElement('att');
att.setAttribute('type','string');
att.setAttribute('name','sbml name');
att.setAttribute('value','Whole_Cell_Model');
graph.appendChild(att);

end

function graph = add_molecule_nodes(docNode, graph, archive, options)

% Total number of molecules in stoichiometry structure
tMolecule = length(archive.stoichiometry.molecule);

% Number of molecule nodes in layout list
try
    ttMolecule = length(options.layout.molecules);
catch
    ttMolecule = 0;
end

% For every molecule in the stoichiometry list
for iMolecule = 1:tMolecule
    
    % Progress
    display_progress(iMolecule,tMolecule);
    
    % Information about molecule
    name = archive.stoichiometry.molecule(iMolecule).name;
    tag = archive.stoichiometry.molecule(iMolecule).tag;
    compartment = archive.stoichiometry.molecule(iMolecule).compartment;
    id = archive.stoichiometry.molecule(iMolecule).id;
    type = archive.stoichiometry.molecule(iMolecule).type;
    number = find(strcmp({archive.info.(type).ID},id));
    try
        category = archive.info.Metabolite(number).Category;
    catch
        category = '';
    end
    
    % Find match inside the layout
    idMatch = arrayfun(@(x)(strcmpi(options.layout.molecules(x).wid,id)),1:ttMolecule);
    compMatch = arrayfun(@(x)(strcmpi(options.layout.molecules(x).c_wid,compartment)),1:ttMolecule);
    match = find(idMatch .* compMatch == 1);
    
    % Total number of nodes to be made
    tNode = length(match);
    
    % If there is no molecule in layout: create a new molecule
    if tNode == 0
        tNode = 1;
        new = true;
    else
        new = false;
    end
    
    % Make a node for every match
    for iNode = 1:tNode
        
        % Node specifi tag
        specific_tag = [tag '_' num2str(iNode)];
        
        % Attributes
        node = docNode.createElement('node');
        node.setAttribute('label',specific_tag);
        node.setAttribute('id',specific_tag);
        node.setAttribute('raw_id',tag);
        
        % Attributes list
        attributeList = {'string','canonicalName',specific_tag;...
            'string','compartment',compartment;...
            'string','id',id;...
            'string','name',name;...
            'string','type',type;...
            };
        
        % Category field if exist
        if isfield(archive.info.Metabolite,'Category')
            attributeList(end+1,:) = {'string','category',category};
        end
        
        % Set attributes
        node = set_attributes(docNode, node, attributeList);
        
        % Get grpahics structure
        graphicsStructure = options.moleculeNode;
        
        % Apply layout to graphics structure if not new
        if ~new
            graphicsStructure.x = options.layout.molecules(match(iNode)).x;
            graphicsStructure.y = options.layout.molecules(match(iNode)).y;
        end
        
        % Set graphics
        node = set_graphics(docNode, node, graphicsStructure, iMolecule);
        
        % Add molecule node to graph node
        graph.appendChild(node);
    end
    
end

end

function graph = add_reaction_nodes(docNode, graph, archive, options)

% Total number of reactions in stoichiometry structure
tReaction = length(archive.stoichiometry.reaction);

% Number of reaction nodes in layout
try
    ttReaction = length(options.layout.reactions);
catch
    ttReaction = 0;
end

% For every reaction in the stoichiometry list
for iReaction = 1:tReaction
    
    % Progress
    display_progress(iReaction,tReaction);
    
    % Information about molecule
    name = archive.stoichiometry.reaction(iReaction).name;
    reversible = archive.stoichiometry.reaction(iReaction).reversible;
    reactionFormula = archive.stoichiometry.reaction(iReaction).reactionFormula;
    id = archive.stoichiometry.reaction(iReaction).id;
    number = find(strcmp({archive.info.MetabolicReaction.ID},id));
    
    % Find match inside the layout
    idMatch = arrayfun(@(x)(strcmpi(options.layout.reactions(x).wid,id)),1:ttReaction);
    match = find(idMatch);
    
    % Total number of nodes to be made
    tNode = length(match);
    
    % If there is no reaction match in layout: create a new reaction node
    if tNode == 0
        tNode = 1;
        new = true;
    else
        new = false;
    end
    
    % Make a node for every match
    for iNode = 1:tNode
        
        % Node specifi tag
        specific_id = [id '_' num2str(iNode)];
        
        % Attributes
        node = docNode.createElement('node');
        node.setAttribute('label',specific_id);
        node.setAttribute('id',specific_id);
        node.setAttribute('raw_id',id);
        
        % Attributes list
        attributeList = {'string','canonicalName',specific_id;...
            'string','id_original',id;...
            'string','reaction name',name;...
            'string','reaction number',number;...
            'string','reaction formula',reactionFormula;...
            'boolean','reversible',reversible;...
            'string','type','MetabolicReaction'};
        
        % Set attributes
        node = set_attributes(docNode, node, attributeList);
        
        % Get grpahics structure
        graphicsStructure = options.reactionNode;
        
        % Apply layout to graphics structure if not new
        if ~new
            graphicsStructure.x = options.layout.reactions(match(iNode)).x;
            graphicsStructure.y = options.layout.reactions(match(iNode)).y;
        end
        
        % Set graphics
        node = set_graphics(docNode, node, graphicsStructure, iReaction);
        
        % Add molecule node to graph node
        graph.appendChild(node);
    end
    
end

end

function graph = add_protein_nodes(docNode, graph, archive, options)

% Total number of protein in stoichiometry structure
tProtein = length(archive.stoichiometry.protein);

% Number of molecule nodes in layout
try
    ttMolecule = length(options.layout.molecules);
catch
    ttMolecule = 0;
end

% For every reaction in the stoichiometry list
for iProtein = 1:tProtein
    
    % Progress
    display_progress(iProtein,tProtein);
    
    % Information about protein
    name = archive.stoichiometry.protein(iProtein).name;
    tag = archive.stoichiometry.protein(iProtein).tag;
    compartment = archive.stoichiometry.protein(iProtein).compartment;
    id = archive.stoichiometry.protein(iProtein).id;
    tmp1 = archive.stoichiometry.protein(iProtein).reaction; % There might be multiple reactions affected
    tmp2 = [arrayfun(@(x)({[tmp1{x} ':']}),1:length(tmp1))]; % Make a list that includes saperators ":"
    reaction = [tmp2{:}]; % Convert to one string
    type = archive.stoichiometry.protein(iProtein).type;
    number = find(strcmp({archive.info.(type).ID},id));
    number = number(1); % Select the first match
    
    % Find match inside the layout
    idMatch = arrayfun(@(x)(strcmpi(options.layout.molecules(x).wid,id)),1:ttMolecule);
    match = find(idMatch);
    
    % Total number of nodes to be made
    tNode = length(match);
    
    % If there is no reaction match in layout: create a new reaction node
    if tNode == 0
        tNode = 1;
        new = true;
    else
        new = false;
    end
    
    % Make a node for every match
    for iNode = 1:tNode
        
        % Node specifi tag
        specific_id = [id '_' num2str(iNode)];
        specific_tag = [tag '_' num2str(iNode)];
        
        % Attributes
        node = docNode.createElement('node');
        node.setAttribute('label',specific_tag);
        node.setAttribute('id',specific_tag);
        node.setAttribute('raw_id',tag);
        
        % Attributes list
        attributeList = {'string','canonicalName',specific_id;...
            'string','id_original',id;...
            'string','tag',tag;...
            'string','compartment',compartment;...
            'string','reaction name',name;...
            'string','protein number',number;...
            'string','reaction',reaction;...
            'string','type',type};
        
        % Set attributes
        node = set_attributes(docNode, node, attributeList);
        
        % Get grpahics structure
        graphicsStructure = options.proteinNode;
        
        % Apply layout to graphics structure if not new
        if ~new
            graphicsStructure.x = options.layout.molecules(match(iNode)).x;
            graphicsStructure.y = options.layout.molecules(match(iNode)).y;
        end
        
        % Set graphics
        node = set_graphics(docNode, node, graphicsStructure, iProtein);
        
        % Add molecule node to graph node
        graph.appendChild(node);
    end
    
end

end

function graph = add_reaction_edges(docNode, graph, archive, options, lNode)

% Total number of reactions in stoichiometry structure
tReaction = length(archive.stoichiometry.reaction);

% For every reaction
for iReaction = 1:tReaction
    
    % Progress
    display_progress(iReaction,tReaction);
    
    % Get id of reaction
    idReaction = archive.stoichiometry.reaction(iReaction).id;
    
    % Find all the reaction nodes with the same id
    lReactionNode = find(strcmp({lNode.raw_id},idReaction));
    
    % Get information about the molecules involved in the reaction
    infoMolecule = archive.stoichiometry.reaction(iReaction).molecule;
    
    % total number of reaction nodes
    tReactionNode = length(lReactionNode);
    
    % For every reaction node: make a edge between the reaction and the
    % closest moleculeNode for every molecule involved in the reaction
    for iReactionNode = 1:tReactionNode
        
        % Node number of reaction
        iNode = lReactionNode(iReactionNode);
        
        % Find the closest node of each molecule that is involved in the
        % reaction
        targetNode = closest_node(lNode, iNode, {infoMolecule.tag});
        
        % Total number of molecule nodes
        tMoleculeNode = length(targetNode);
        
        for iMoleculeNode = 1:tMoleculeNode
            
            % Get the participation
            participation = archive.stoichiometry.reaction(iReaction).molecule(iMoleculeNode).participation;
            
            % Determine source and target
            switch participation
                case 'reactant'
                    source = lNode(targetNode(iMoleculeNode)).id;
                    raw_source = lNode(targetNode(iMoleculeNode)).raw_id;
                    target = lNode(iNode).id;
                    raw_target = lNode(iNode).raw_id;
                    labelElement = '(reactant-reaction)';
                case 'product'
                    source = lNode(iNode).id;
                    raw_source = lNode(iNode).raw_id;
                    target = lNode(targetNode(iMoleculeNode)).id;
                    raw_target = lNode(targetNode(iMoleculeNode)).raw_id;
                    labelElement = '(reaction-product)';
            end
            
            % Make label
            label = [source labelElement target];
            
            % Create edge
            edge = docNode.createElement('edge');
            edge.setAttribute('label',label);
            edge.setAttribute('source',source);
            edge.setAttribute('target',target);
            edge.setAttribute('raw_source',raw_source);
            edge.setAttribute('raw_target',raw_target);
            graph.appendChild(edge);
            
            % Attributes list
            attributeList = {'string','canonicalName',label;...
                'string','interaction',labelElement;...
                'real','sbml stoichiometry','1.0'};
            
            % Set attributes
            edge = set_attributes(docNode, edge, attributeList);
            
            % Graphics
            edge = set_graphics(docNode, edge, options.reactionEdges, iNode);
            
        end
        
    end
    
end

end

function graph = add_protein_edges(docNode, graph, archive, options, lNode)

% Total number of proteins in stoichiometry structure
tProtein = length(archive.stoichiometry.protein);

% For every protein
for iProtein = 1:tProtein
    
    % Progress
    display_progress(iProtein,tProtein);
    
    % Get id of protein
    idProtein = archive.stoichiometry.protein(iProtein).tag;
    
    % Find all the protein nodes with the same id
    lProteinNode = find(strcmp({lNode.raw_id},idProtein));
    
    % Get information about the reactions affected by the protein
    infoProtein = archive.stoichiometry.protein(iProtein);
    
    %     % Total number of protein nodes
    %     tProteinNode = length(lProteinNode);
    
    % Total number of reaction types affected by protein
    tReaction = length(infoProtein.reaction);
    
    % Check for all the reaction types affected by the protein
    for iReaction = 1:tReaction
        
        % Get the idea of the reaction affected by the protein
        idReaction = infoProtein.reaction{iReaction};
        
        % Find all the reaction nodes with the same id
        lReactionNode = find(strcmp({lNode.raw_id},idReaction));
        
        % Total number of reaction nodes
        tReactionNode = length(lReactionNode);
        
        % For every reaction node: make a edge between the reaction and the
        % closest proteinNode
        for iReactionNode = 1:tReactionNode
            
            % Node number of reaction
            iNode = lReactionNode(iReactionNode);
            
            % Find the closest node of each reaction that is affected by the
            % protein
            targetNode = closest_node(lNode, iNode, {idProtein});
            
            % Total number of molecule nodes
            tMoleculeNode = length(targetNode);
            
            % Source and target info of edge
            source = lNode(targetNode).id;
            raw_source = lNode(targetNode).raw_id;
            target = lNode(iNode).id;
            raw_target = lNode(iNode).raw_id;
            labelElement = '(protein-reaction)';
            
            % Make label
            label = [source labelElement target];
            
            % Create edge
            edge = docNode.createElement('edge');
            edge.setAttribute('label',label);
            edge.setAttribute('source',source);
            edge.setAttribute('target',target);
            edge.setAttribute('raw_source',raw_source);
            edge.setAttribute('raw_target',raw_target);
            graph.appendChild(edge);
            
            % Attributes list
            attributeList = {'string','canonicalName',label;...
                'string','interaction',labelElement;...
                'real','sbml stoichiometry','1.0'};
            
            % Set attributes
            edge = set_attributes(docNode, edge, attributeList);
            
            % Graphics
            edge = set_graphics(docNode, edge, options.proteinEdges, iNode);
            
        end
        
    end
    
end

end

function newLayout = format_layout(layout)
% Change format of layout structure

% Get fields
lField = fields(layout);

% Number of fields
tField = length(lField);

for iField = 1:tField
    
    % Create field
    newLayout.(lField{iField}) = {};
    
    % Number of elements
    tElement = length(layout.(lField{iField}));
    
    for iElement = 1:tElement
        if iElement ~= 1
            newLayout.(lField{iField})(end+1) = layout.(lField{iField}){iElement};
        else
            newLayout.(lField{iField}) = layout.(lField{iField}){iElement};
        end
        
    end
    
end

end

function node = set_attributes(docNode, node, attributeList)
% Set attributes from attributes list to node
% Note: attribute with name 'interaction' will be handled differently

% Total number of attributes
tAttribute = length(attributeList);

% Add attributes
for iAttribute = 1:tAttribute
    
    % Value
    type = attributeList(iAttribute,1);
    name = attributeList(iAttribute,2);
    value = attributeList(iAttribute,3);
    
    % Check if value matches type
    if ~ischar(value{:}) && strcmp(type{:},'string')
        value = {'error: no value found'};
    elseif ~ischar(value{:}) && strcmp(type{:},'real')
        value = {num2str(value{:})};
    end
    
    % Add attribute to node
    att = docNode.createElement('att');
    att.setAttribute('type',type);
    att.setAttribute('name',name);
    att.setAttribute('value',value);
    % Exception:
    if strcmp(attributeList(iAttribute,2),'interaction')
        att.setAttribute('cy:editable','false');
    end
    node.appendChild(att);
end
end

function node = set_graphics(docNode, node, structure, iNode)
% Add graphics of options to node

% Get graphics fields
lField = fields(structure);

% Initiate graphics
graphics = docNode.createElement('graphics');

% Total number of graphics
tGraphic = length(lField);

% Add graphics
for iGraphic = 1:tGraphic
    
    % Graphics field
    graphic = lField{iGraphic};
    
    % Graphics value
    value = structure.(graphic);
    
    % Change format of value
    if ~ischar(value)
        value = num2str(value);
    elseif strcmp(value(1:6),'shift_')
        value = str2double(value(7:end)) * iNode;
        value = num2str(value);
    end
    
    % Add value and graphic to node
    graphics.setAttribute(graphic,value);
    node.appendChild(graphics);
end
end

function lNode = process_nodes(graph)
% Get location, id, raw_id, number information of all the existing nodes

% Get all the current nodes (id, raw_id and number)
lNode = getTypeOfTarget(graph,'node',{'id','raw_id'});

% Total number of nodes
tNode = length(lNode);

% Process every node
for iNode = 1:tNode
    
    % Open node
    node = item(graph,lNode(iNode).number);
    
    % x and y values of graphics
    lGraphic = getTypeOfTarget(node,'graphics',{'x','y'});
    
    % Add x and y values to node list
    lNode(iNode).x = lGraphic.x;
    lNode(iNode).y = lGraphic.y;
    
end

end

function targetNode = closest_node(lNode, iNode, molecules)
% targetNode = closest_node(lNode, nodeNumber, targets)
%
%   lNode is the list of all the nodes which includes the x, y, id, raw_id,
%   and number of the nodes
%   iNode is the node number which should be linked
%   molecules is the list of tags/id's that should be linked with the node

% Total number of molecules involved in reaction
tMolecule = length(molecules);

% Make an edge between every reaction and one of each molecule type
% Node
for iMolecule = 1:tMolecule
    
    % Find all the reaction nodes with the same id
    lMoleculeNode = find(strcmp({lNode.raw_id},molecules{iMolecule}));
    
    % Find closest Molecule node of every type of molecule
    
    % Find distance between reaction node and molecule nodes
    diffX = str2double(lNode(iNode).x) - str2double({lNode(lMoleculeNode).x});
    diffY = str2double(lNode(iNode).y) - str2double({lNode(lMoleculeNode).y});
    
    % Find closest node
    [~,closestNode] =  min(diffX.^2 + diffY.^2);
    targetNode(iMolecule) = lMoleculeNode(closestNode);
    
end

end