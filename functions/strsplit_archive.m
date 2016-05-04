% An extended version of strsplit which works the same as both strsplit of
% the Matlab toolbox and whole cell model

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/02/2016

function varargout = strsplit_archive(varargin)

% Investigate function
type = which('strsplit');

% Check if it's the strfind of the Matlab Toolbox or Whole cell model
% Empty: Matlab toolbox; Not empty: Whole cell model
check = strfind(type,'lib/util/strutil');

if isempty(check)
    % Matlab version
    out = strsplit(varargin{:});
else
    % Whole cell model version
    if nargin > 1
        out = strsplit(varargin{2},varargin{1});
    else
        out = strsplit(' ',varargin{1});
    end
    
end

% Set output
varargout{1} = out;

end