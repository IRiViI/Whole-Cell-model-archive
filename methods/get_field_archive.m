% Get structure of parent structure
%   Example:
%   archive.info.x.y, where length(archive.info.x) > 1 and you want all the
%   fields y of structure archive.info.x.
%   y = get_field_archive(archive, {'info','x'}); where length(y) =
%   length(archive.info.x)
%   
%   Alternative
%   you want the the same of the second info structure archive.info(2).x.y.
%   In that case you have to specify the route:
%   y = get_field_archive(archive,{'info','x'},[2,0]);
%   Where 2 stands for the second info structure and 0 for all the x
%   structures

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/02/2016

function field = get_field_archive(archive, varargin)

% Specifiy the route
field = archive;
nField = varargin{1};
if nargin > 2
    lNumber = varargin{2};
else
    lNumber = zeros(size(nField));
end

% Get to field
for iField = 1:length(nField)
    field = field.(nField{iField});
    if lNumber(iField) == 0
    else % Select specific field
        field = field(lNumber(iField));
    end
    read_length(iField) = length(field);
end

end