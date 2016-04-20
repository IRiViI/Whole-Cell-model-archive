function archive = add_sets_archive(archive, varargin)
%This functions adds sets to the archive object. It looks for the
%   options.mat files within the folder and its subfolders.
%   The information about the folders and seeds of these simulations are
%   ordered accordingly in the output structure.
%
%   archive = archive_list(archive, 'path/to/folder')
%   This directory should contain the "options.mat" file of the simulation
%   in order to execute this function. The other files that should be
%   analysed using the output list should be in the same folder as the
%   "options.mat" that belongs to that specific simulation. These files
%   could still be added to their folder at a later point of time.
%
%   'progress'      - Show progress of function. ('On' or 'Off'). Default
%                     setting is 'On'.
%  'knockout_filter'- Make a new set when the knockouts do not matches the
%                     knockouts of a previous set. Default is 'On'

% Mendatory inputs
if nargin > 0
    file = varargin{1};
else
    error('Not enough input arguments. Please give the path to the folder containing the simulation files');
end

% Default options
options.progress = 'On';
options.knockout_filter = 'On';
options.short = 'Off';

% Adjust default options
if nargin > 2
    input_options = struct(varargin{2:end});
    field_names = fieldnames(input_options);
    for i = 1:size(field_names,1)
        options.(field_names{i}) = input_options.(field_names{i});
    end
end

% options for find_files function ('short' will result in a less thorough
% search in the directory (much faster).
ffstruct = {};
if strcmp(options.short,'On')
    ffstruct{end+1} = 'short';
end

% If there is no propperty "set" yet (Old dated as well, because it is always present)
if ~isprop(archive,'set')
    archive.addprop('set');
end

% Find all the option files
file_path_structure = find_files(file,'options.mat','On',ffstruct{:});

% Progress spacer
if strcmp(options.progress,'On')
    fprintf('Progress set simulations to archive:      ');
end

% Number of simulation files
t_simulation = length(file_path_structure);

% Keep track of how many fields are missing
missing_fields_sim = {};
missing_fields_set = {};

% For every simulation file
for i_simulation = 1:t_simulation
    
    % Show progress
    if strcmp(options.progress,'On')
        display_progress(i_simulation,t_simulation)
    end
    
    % Path to options file
    options_file_path = [file_path_structure(i_simulation).folder '/' file_path_structure(i_simulation).file];
    
    % Load the options structure of the simulation.
    options_structure = load(options_file_path,'seed','geneticKnockouts');
    
    % Get knockout information
    knockout = options_structure.geneticKnockouts;
    if isempty(knockout)
        knockout = {''};       % Name, if no knockouts is found
    end
    t_knockout = length(knockout);
    
    % Get the total number of sets
    t_set = length(archive.set);
    % If nothing happens, the number for the set is one greater than
    % the current newest set.
    set_number = t_set + 1;
    
    % Check all the knockouts if you want to make different sets for
    % different gene disrupted strains
    if strcmp(options.knockout_filter,'On')
        
        % Check if there are already some sets present in the archive list
        if isfield(archive.set,'knockout')
            % If so, make a list with all the knockouts of every set
            knockouts = archive.extract_data({'knockout'},'select','all');
        else
            knockouts = [];
        end
        % Check if the knockout combination already exists
        for i_set = 1:t_set
            % Get the number of knockouts for exist set "i_set"
            t_knockout_i_set = length(knockouts{i_set});
            % Continue if the number of knockouts match
            if t_knockout_i_set == t_knockout
                % Make a check list with the same length as the number of
                % knockouts
                check = zeros(1,t_knockout_i_set);
                % Check for every gene knockout of the simulation if it also
                % knocked out in set "i_set" and use the "check" quantity to
                % keep track of the gene knockouts
                for i_knockout = 1:t_knockout
                    % Check knockout "i_knockout" and see if it's present in
                    % set "i_set"
                    check = check + strcmp(knockouts{i_set},knockout(i_knockout));
                end
                % If A perfectmatch is found, change the set_number to the
                % i_set number
                if sum(check) == t_knockout && sum(check) == t_knockout_i_set && max(check) == 1
                    set_number = i_set;
                    break;
                end
            end
        end
    end
    
    % Make new simulation structure
    simulation = struct();
    simulation.seed = options_structure.seed;
    simulation.folder = file_path_structure(i_simulation).folder;
    
    % Check if this is a new type of set or if it should be placed in an
    % preexisting set
    if set_number == t_set+1 || isempty(archive.set)
        
        % Create a new set
        set = struct();
        set.knockout = knockout;
        set.simulation = simulation;
        if isempty(archive.set)
            % If this will be the first set of the archive structure:
            % Add set:
            archive.set = set;
        else
            % If there are already sets present in the archive structure:
            % Make set compatable with the preexisting sets:
            field = fields(archive.set(1));
            tField = length(field);
            for iField = 1:tField
               if ~isfield(set,field(iField))
                   set.(field{iField}) = [];
                   missing_fields_set(end+1) = field(iField);
               end
            end
            % Add set:
            archive.set(set_number) = set;
        end
    else
        % If there are already simulations present in the set structure:
            % Make set compatable with the preexisting simulations:
            field = fields(archive.set(set_number).simulation(1));
            tField = length(field);
            for iField = 1:tField
               if ~isfield(simulation,field(iField))
                   simulation.(field{iField}) = [];
                   if sum(strcmp(missing_fields_sim,field(iField)))<1 ||...
                           isempty(missing_fields_sim)
                       missing_fields_sim(end+1) = field(iField);
                   end
               end
            end
        % Add simulation to existing set
        archive.set(set_number).simulation(end+1) = simulation;
    end
    
end

% Update if there are some fields missing
if length(missing_fields_set) > 0
    fprintf('The new set(s) are missing %d field(s)\n',length(missing_fields_set))
end
if length(missing_fields_sim) > 0
    fprintf('The new sim(s) are missing %d field(s)\n',length(missing_fields_sim))
end

end