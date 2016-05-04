% Return whether the string ends with ".mat" or not. 
%
% state = check_direct('file.mat')
%   state = true or false

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/02/2016

function state = check_direct(string)

% Separator
symbol = '.';

% Get elements string
C = strsplit_archive(string,symbol);

% Check the elements after the "."
if strcmp(C(end),'mat')
    state = true;
else
    state = false;
end

end