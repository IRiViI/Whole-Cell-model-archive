
% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/02/2016

function field = get_field_archive(archive, varargin)

% Get initial structure
field = archive;
nField = varargin{1};
if nargin > 1
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