% Convert the layout generated by cytoscape (export json file (Network and 
%   View)) to the format used in function "cytoscape_xml_list.m".
%
%   json_layout_converter('path/to/json_file.cyjs')
%
%   Optional:
%   'out'       - Set the path or name of new json file
%   'json_lib'  - Set the path to the json library
%
%   NOTE: This function uses the "jsonlab" library. Url:
%   http://www.mathworks.com/matlabcentral/fileexchange/
%   33381-jsonlab--a-toolbox-to-encode-decode-json-files

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/13/2016

function json_layout_converter(varargin)

% Mandatory inputs
options.in = varargin{1};

% Default settings
options.json_lib = 'functions/others/jsonlab';

% Change default settings
if nargin > 1
    input_options = struct(varargin(2:end));
    field_names = fieldnames(input_options);
    for i = 1:size(field_names,1)
        options.(field_names{i}) = input_options.(field_names{i});
    end
end

% ------------- Process request ------------- %

% Output name
if ~isfield(options,'out')
    tmp = strsplit(options.in,'/');
    options.out = ['output' '/' 'converted_' tmp{end}];
end

% Add path to json library
addpath(options.json_lib);

% Get layout
layout_in = loadjson(options.in);

% Initate out layout
layout_out = struct();

% Total number of nodes
t_node = length(layout_in.elements.nodes);

% Number of reactions and molecules
i_react = 0;
i_molec = 0;

% Progress spacing
fprintf('Progress: \n     ');

% For every node
for i_node = 1:t_node
    
    % Show progress
    display_progress(i_node,t_node)
    
    % Get node
    node = layout_in.elements.nodes{i_node};
    
    % If molecule/metabolite
    if strcmpi(node.data.type,'Metabolite') ||...
            strcmpi(node.data.type,'ProteinComplex') ||...
            strcmpi(node.data.type,'ProteinMonomer') ||...
            strcmpi(node.data.type,'Stimulus')
        
        i_molec =  i_molec + 1;
        
        layout_out.molecules{i_molec}.wid = ...
            node.data.id_original;
        if isempty(node.data.id_original) % Once, something went wrong and I have(had) to obtain the id from the name...
            tmp = strsplit(node.data.name,'_');
            tmp2 = [arrayfun(@(x)([tmp{x} '_']),1:length(tmp)-1,'UniformOutput',false)];
            tmp3 = [tmp2{:}];
            id = tmp3(1:end-1);
            layout_out.molecules{i_molec}.wid = id;
        end
        
        layout_out.molecules{i_molec}.c_wid = ...
            node.data.compartment;
        
        layout_out.molecules{i_molec}.x = ...
            node.position.x;
        
        layout_out.molecules{i_molec}.y = ...
            node.position.y;
        
    % If reaction
    elseif strcmpi(node.data.type,'MetabolicReaction')
        i_react =  i_react + 1;
        
        layout_out.reactions{i_react}.wid = ...
            node.data.id_original;
        if isempty(node.data.id_original) % If the original id does not exist, see if you can get the id from the name
            tmp = strsplit(node.data.name,'_');
            tmp2 = [arrayfun(@(x)([tmp{x} '_']),1:length(tmp)-1,'UniformOutput',false)];
            tmp3 = [tmp2{:}];
            id = tmp3(1:end-1);
            layout_out.reactions{i_react}.wid = id;
        end
        
        layout_out.reactions{i_react}.x = ...
            node.position.x;
        
        layout_out.reactions{i_react}.y = ...
            node.position.y;
    else
        warning_message = sprintf('feature type "%s" of the feature %s is not recognized',node.data.type,node.data.id_original);
        warning(warning_message);
    end
end

% Save json file
savejson('',layout_out,options.out);

% Progresss
fprintf(['Done: ' options.out '\n'])

end