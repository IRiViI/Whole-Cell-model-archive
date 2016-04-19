% Extract a certain value out of the archive. 
% value = extract_data_archive(archive,{'cluster','value'})
%
%   'set'       - The set number (only required when the target is
%                 simulation)
%   'select'    - The coordinates of the value to extract. If the value
%                 is within in a vector or matrix, specify the x and y 
%                 values Example: list.data(1).field.subfield = rand(3,3)).
%                 'select',[2,2],... selects the element in the middle. 
%                 Other options are "'select','all'" and "'select','end'".

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 11/03/2016

function value = extract_data_archive(varargin)

% Mandatory input
archive = varargin{1};
path = varargin{2}; % Example: path = {'cluster','values'};

% Change the format if the format is "'field'" instead of "{'field'}"
if isstr(path)
   path = {path};
end

% Default settings
% options.select = ones(size(path));

% Target the strains or the simulation of a strain
options.target = 'strain'; % 'strain', 'simulation'
options.set = 1;
options.select = 'all';

% Optional inputs
if nargin > 2
    % Get inputs
    input_options = struct(varargin{3:end});
    
    % Number and names of input options
    nOptions = fields(input_options);
    tOptions = length(nOptions);
    
    % Apply option adjustments
    for iOptions = 1:tOptions
        options.(nOptions{iOptions}) = input_options.(nOptions{iOptions});
        if strcmp(nOptions{iOptions},'set')
            options.target = 'simulation';
        end
    end
    
end

value = [];

% Assign select
% select = options.select;

% Get values of a strain data request
if strcmp(options.target,'strain')
    
    % Total number of strains
    tStrain = length(archive.set);
    
    % For every strain
    for iStrain = 1:tStrain
        
        % If there is a value, get value
        try
            switch length(path)
                case 1
                    out = archive.set(iStrain)...
                        .(path{1});
                case 2
                    out = archive.set(iStrain)...
                        .(path{1})...
                        .(path{2});
                case 3
                    out = archive.set(iStrain)...
                        .(path{1})...
                        .(path{2})...
                        .(path{3});
                case 4
                    out = archive.set(iStrain)...
                        .(path{1})...
                        .(path{2})...
                        .(path{3})...
                        .(path{4});
            end
            
            % Select
            if strcmp(options.select,'all')
            elseif strcmp(options.select,'end')
                out = out(end);
            else
                out = out(options.select);
            end
            
            % Save value to value array
            if isnumeric(out)
                if length(out) < 2
                    value(iStrain) = out;
                else
                    value{iStrain} = out;
                end
            else
                value{iStrain} = out;
            end
            
            % If there is no value found
        catch
            value(iStrain) = nan;
        end
        
    end
    
elseif strcmp(options.target,'simulation')
    
    % Total number of simulations
    tSimulation = length(archive.set(options.set).simulation);
    
    % For every simulation
    for iSimulation = 1:tSimulation
        
        % If there is a value, get value
        try
            switch length(path)
                case 1
                    out = archive.set(options.set).simulation(iSimulation)...
                        .(path{1});
                case 2
                    out = archive.set(options.set).simulation(iSimulation)...
                        .(path{1})...
                        .(path{2});
                case 3
                    out = archive.set(options.set).simulation(iSimulation)...
                        .(path{1})...
                        .(path{2})...
                        .(path{3});
                case 4
                    out = archive.set(options.set).simulation(iSimulation)...
                        .(path{1})...
                        .(path{2})...
                        .(path{3})...
                        .(path{4});
            end
            
            % Select
            if strcmp(options.select,'all')
            elseif strcmp(options.select,'end')
                out = out(end);
            else
                out = out(options.select);
            end
            
            % Save value to value array
            if isnumeric(out)
                if length(out) < 2
                	value(iSimulation) = out;
                else
                    value{iSimulation} = out;
                end
            else
                value{iSimulation} = out;
            end
            
            % If there is no value found
        catch
            value(iSimulation) = nan;
        end
        
    end
    
end

end

