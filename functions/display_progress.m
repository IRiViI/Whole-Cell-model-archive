%Display progress of function in command window.
%   display_progress(current_number,total_number)
%
%   Example:
%   display_progress(1,10)
%       Display: 10%
%
%   Hint: Start by adding 5 spaces to the command window before running
%   this script for the first time.


% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 01/29/2016

function display_progress(varargin)

% Process mendatory input arguments
cN = varargin{1};
tN = varargin{2};

% Calculate the percentage done
pN = ceil(cN/tN*100);

% Add spacer (to correct for string length)
if pN >= 100
    spacer =  '';
elseif pN >= 10
    spacer =  ' ';
else
    spacer =  '  ';
end

% Update progress
fprintf('\b\b\b\b\b%d%%%s\n',pN,spacer);

end
