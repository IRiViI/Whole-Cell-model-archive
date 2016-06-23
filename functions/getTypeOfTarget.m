%Find the elements with a specific attribute in a xml structure.
%   info = getTypeOfTarget(doc,target,type)
%   The type can be both a string or a char string
%   Examples:
%       info = getTypeOfTarget(doc,'node','id')
%       info = getTypeOfTarget(doc,'node',{'id','name'})

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 02/05/2016

function info = getTypeOfTarget(varargin)

% Proces input arguments
doc = varargin{1};
target = varargin{2};
type = varargin{3};

% Get the number of "att" and "node" elements
dNodes = getChildNodes(doc);
tNode = getLength(dNodes);

% Create info list
info = struct();

% Entry counter
iEntry = 1;

for iNode = 0:tNode-1
    
    % Get node
    node = item(doc,iNode);
    
    % Name of the nod
    name = toCharArray(getNodeName(node))';
    
    % If it's the
    if strcmp(name,target)
        
        % Get the value(s)
        if iscellstr(type)  % If it's a char string
            tType = length(type);
            for iType = 1:tType
                value = toCharArray(getAttribute(node,type{iType}))';
                info(iEntry).(type{iType}) = value;
            end
        else    % If it's a string
            value = toCharArray(getAttribute(node,type))';
            info(iEntry).(type) = value;
        end
        
        % Get the number of this node
        info(iEntry).number = iNode;
        
        % Update entry counter
        iEntry = iEntry + 1;
        
    end
    
end

end