% T-test between two sets of sets/simulations (see ttest2 for more details 
%   about ttest)
%   [H,P,CI] = snap_shot_Ttest_archive(archive, snapShot, set1, set2)
%   snapShot is snap shot number of archive structure
%   set1 is an array of sets for the first set t-test
%   set2 is an array of sets for the second set of t-test
%
%   'simulation1'   - Select a smaller set of the simulation to include for
%                     set 1
%   'simulation2'   - Select a smaller set of the simulation to include for
%                     set 2
%   'ttvarargin'    - varargin settings for ttest2 function Enclose options
%                     in withbetween double curly brackets "{{...}}". See
%                     example
%
%   Examples:
%       [H,P] = snap_shot_Ttest_archive(archive,3,1,[2,3,4])
%               % Compare all the simulations of sets 2,3 and 4 as one set
%               with set 1. 
%
%       [H,P] = snap_shot_Ttest_archive(archive,3,1,2,...
%               'ttvarargin',{{'alpha',0.05}})

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/25/2016

function varargout = snap_shot_Ttest_archive(archive, snapShot, set1, set2, varargin)
%% Settings

% Default settings
options.set1 = set1;
options.set2 = set2;
options.snapShot = snapShot;
options.simulation1 = 0; 	% Default: analyze all simulations
options.simulation2 = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.ttvarargin = {};

% Change settings
tmpOptions = struct(varargin{1:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

%% Process settings

%% Run code

% Create the first matrix
matrix1 = snap_shot_matrix(archive, options.snapShot,options.set1,options.simulation1);

% Create the second matrix
matrix2 = snap_shot_matrix(archive, options.snapShot,options.set2,options.simulation2);

tElement = size(matrix1,1);
tCompartment = size(matrix1,2);
a = []; b = []; c = [];

for iElement = 1:tElement
    
    for iCompartment = 1:tCompartment
        
        [a(iElement,iCompartment,:),...
            b(iElement,iCompartment,:),...
            c(iElement,iCompartment,:),...
            d(iElement,iCompartment,:)] =...
            ttest2(squeeze(matrix1(iElement,iCompartment,:)),...
            squeeze(matrix2(iElement,iCompartment,:)),...
            options.ttvarargin{:});
        
    end
    
end

% Display progress
if options.progress
    fprintf('Done\n');
end

%% Set outputs

% Set output
varargout{1} = a;
varargout{2} = b;
varargout{3} = c;
varargout{4} = d;

end

function matrix = snap_shot_matrix(archive,snapShot,lSet,lSimulation)

% Initiate matrix structure
structure = struct([]);

% Process sets
for iSet = lSet
    
    % Make list of sets to analyze
    if lSimulation == 0
        tSimulation = archive.simulations(iSet);    % Number of simulations to analyze
        lSimulation = 1:tSimulation;                % List of simulations to analyze
    else
        lSimulation = lSimulation;                  % List of simulations to analyze
        tSimulation = length(lSimulation);          % Number of simulations to analyze
    end
    
    
    % Process simulations
    for iSimulation = lSimulation
        
        structure(end+1).values = archive.set(iSet).simulation(iSimulation).snapShot(snapShot).values;
        
    end
    
end

% Create combined matrix
matrix = cat(3 + 1,structure.values);

end
