function archive = add_sets_archive(archive, varargin)
%This functions adds sets to the archive object. It looks for the
%   options.mat files within the folder and its subfolders.
%   The information about the folders and seeds of these simulations are
%   ordered accordingly in the output structure.
%
%   archive = archive_list(archive, 'path/to/folder')
%   This directory should contain the "options.mat" and a state.mat
%   (like: state-0.mat or mean_compressed_state.mat) file of the simulation
%   in order to execute this function. The other files that should be
%   analysed using the output list should be in the same folder as the
%   "options.mat" that belongs to that specific simulation. These files
%   could still be added to their folder at a later point of time.
%
%   'type'          - label the simulation with a type. Simulations are
%                     only added to a set if their types match.
%   'progress'      - Show progress of function. (true or false). Default
%                     setting is true.
%  'knockout_filter'- Filter sets for matching genetic knockouts
%   'model_filter'  - Filter sets for matching number of
%                     molecules/reactions
%   'new'           - Set 'new' to true if the preexisting sets should be
%                     kept untouched by the newly added simulations
%   'note'          - Add a commend to a set.
%   'short'         - Seriously reduce the time find folder process by
%                     setting 'short' to true. However, it only works if
%                     the simulation folders are at the same level in the
%                     parent folder. Example: Same level: parent/child/sim1 and
%                     parent/child/sim2; Not same level: parent/sim1 and
%                     parent/child/sim2. Default: false

% Mendatory inputs
if nargin > 0
    file = varargin{1};
else
    error('Not enough input arguments. Please give the path to the folder containing the simulation files');
end

% This option is used when the old states whould be kept untouched
options.initial_nSets = length(archive.set);

% Default options
options.progress = true;
options.knockout_filter = true;
options.model_filter = true;
options.type_filter = true;
options.short = false;
options.new = false;
options.note = '';
options.type = '';

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
if options.short
    ffstruct{end+1} = 'short';
end

% If there is no propperty "set" yet (Old dated as well, because it is always present)
if ~isprop(archive,'set')
    archive.addprop('set');
end

% Find all the option files
file_path_structure = find_files(file,'options.mat','On',ffstruct{:});

% Progress spacer
if options.progress
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
    if options.progress
        display_progress(i_simulation,t_simulation)
    end
    
    % Folder dir
    folder = file_path_structure(i_simulation).folder;
    
    % ---- Option file
    
    % Path to options file
    options_file_path = [folder '/' file_path_structure(i_simulation).file];
    
    % Load the options structure of the simulation.
    options_structure = load(options_file_path,'seed','geneticKnockouts');
    
    % ---- State file
    
    % Load a state file
    try
        try
            state_structure = load([folder '/' 'state-0.mat']);
        catch
            try
                state_structure = load([folder '/' 'mean_compressed_state.mat']);
            catch
                state_structure = load([folder '/' 'point_compressed_state.mat']);
            end
        end
    catch
        warningMessage = sprintf('No state files recognized in folder %s',folder);
        warning(warningMessage)
    end
    
    % Get state file information
    info = struct();
    info.ProteinComplex.number = size(state_structure.ProteinComplex.counts,1);
    info.ProteinMonomer.number = size(state_structure.ProteinMonomer.counts,1);
    info.Rna.number = size(state_structure.Rna.counts,1);
    info.Metabolite.number = size(state_structure.Metabolite.counts,1);
    info.MetabolicReaction.number = size(state_structure.MetabolicReaction.fluxs,1);
    
    % Get knockout information
    options.knockout = options_structure.geneticKnockouts;
    if isempty(options.knockout)
        options.knockout = {''};       % Name, if no knockouts is found
    end
    options.t_knockout = length(options.knockout);
    
    % Get the total number of sets
    options.t_set = length(archive.set);
    % If nothing happens, the number for the set is one greater than
    % the current newest set.
    options.set_number = options.t_set + 1;
    
    % Find all the matches during the filetering
    set_match = struct;
    
    % ----------------------  Filters ---------------------- 
    
    % Compare Gene knockouts
    set_match = KO_filter(archive, set_match, options);
    
    % Compare the number of molecules and reactions
    set_match = molecule_reaction_filter(archive, set_match, info, options);
    
    % Compare the type of set (manually set)
    set_match = type_filter(archive, set_match, options);
    
    % --------- Finalizing filtering
    
    % Update set number
    % Check if preexisting simulations should be ignored
    if options.new
        minSet = options.initial_nSets;
    else
        minSet = 1;
    end
    
    % NOTE: tSet is already determined in previous process, include new
    % tSet calculation if this is not true anymore
    
    % Get information of all the filters
    fMatch = fields(set_match);
    tMatch = length(fMatch);
    
    % tSet
    tSet = archive.sets;
    
    % sets to check
    lSet = minSet:tSet;
    
    % Check every targeted set if they match. In other words: Check the
    % data in the different filter sets (set_match.*) and check for every
    % set if its number is present
    for iSet = lSet
        % Try to find a match by all different fileters
        try
            tmp = arrayfun(@(x)(find(set_match.(fMatch{x})==iSet)),1:tMatch);
            % All filters match
            options.set_number = set_match.(fMatch{1})(tmp(1)); % Update the set number
            break % There should only be one match
        catch
            % Not all filters match
        end
    end
    
    
    % ----------------------  Process simulation ---------------------- %
    
    % Make new simulation structure
    simulation = struct();
    simulation.seed = options_structure.seed;
    simulation.folder = file_path_structure(i_simulation).folder;
    
    % Check if this is a new type of set or if it should be placed in an
    % preexisting set
    if options.set_number == options.t_set+1 || isempty(archive.set)
        
        % Create a new set
        set = struct();
        set.note = options.note;
        set.knockout = options.knockout;
        set.simulation = simulation;
        set.info = info;
        set.type = options.type;
        
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
                    shuffled_set.(field{iField}) = [];
                    % Update structing that keeps track of missing fields
                    missing_fields_set(end+1) = field(iField);
                else
                    % The set has be shuffled to match the existing order.
                    % Or else shit goes wrong apparently...
                    shuffled_set.(field{iField}) = set.(field{iField});
                end
            end
            % Add set:
            archive.set(options.set_number) = shuffled_set;
        end
    else
        % If there are already simulations present in the set structure:
        % Make set compatable with the preexisting simulations:
        field = fields(archive.set(options.set_number).simulation(1));
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
        archive.set(options.set_number).simulation(end+1) = simulation;
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

function set_match = molecule_reaction_filter(archive, set_match, info, options)
% Check the number of molecules and reactions

% Compare state information with preexisting sets
if options.model_filter
    % initiate set_match.system
    set_match.system = [];
    % Total number of sets
    tSet = length(archive.set);
    % For every set
    for iSet = 1:tSet
        % Info structure of set
        sInfo = archive.set(iSet).info;
        % Check if the number of molecules and reactions are the same
        if sInfo.ProteinComplex.number == info.ProteinComplex.number &&...
                sInfo.ProteinMonomer.number == info.ProteinMonomer.number &&...
                sInfo.Rna.number == info.Rna.number &&...
                sInfo.Metabolite.number == info.Metabolite.number &&...
                sInfo.MetabolicReaction.number == info.MetabolicReaction.number
            
            % They match
            set_match.system(end+1) = iSet;
        else
            % They do not match
        end
    end
end

end

function set_match = KO_filter(archive, set_match, options)

% Check all the knockouts if you want to make different sets for
% different gene disrupted strains
if options.knockout_filter
    set_match.KO = [];
    
    % Check if there are already some sets present in the archive list
    if isfield(archive.set,'knockout')
        % If so, make a list with all the knockouts of every set
        knockouts = archive.extract_data({'knockout'},'select','all');
        % remove the information of the preexisting sets when they
        % should not be checked
        if options.new
            knockouts(1:options.initial_nSets) = {nan};
        end
    else
        knockouts = [];
    end
    
    % Check if the knockout combination already exists
    for i_set = 1:options.t_set
        % Get the number of knockouts for exist set "i_set"
        t_knockout_i_set = length(knockouts{i_set});
        % Continue if the number of knockouts match
        if t_knockout_i_set == options.t_knockout
            % Make a check list with the same length as the number of
            % knockouts
            check = zeros(1,t_knockout_i_set);
            % Check for every gene knockout of the simulation if it also
            % knocked out in set "i_set" and use the "check" quantity to
            % keep track of the gene knockouts
            for i_knockout = 1:options.t_knockout
                % Check knockout "i_knockout" and see if it's present in
                % set "i_set"
                check = check + strcmp(knockouts{i_set},options.knockout(i_knockout));
            end
            % If A perfectmatch is found, change the options.set_number to the
            % i_set number
            if sum(check) == options.t_knockout && sum(check) == t_knockout_i_set && max(check) == 1
                set_match.KO(end+1) = i_set;
            end
        end
    end
end

end

function set_match = type_filter(archive, set_match, options)
% Check the number of molecules and reactions

% Compare state information with preexisting sets
if options.type_filter
    % initiate set_match.system
    set_match.type = [];
    % Total number of sets
    tSet = length(archive.set);
    % For every set
    for iSet = 1:tSet
        % Get type of set
        if isfield(archive.set(iSet),'type')
            sType = archive.set(iSet).type;
        else
            sType = ''; % If field does not exist, set type empty string
        end
        if ~ischar(sType)
            sType = ''; % If sType is ill defined, set type to empty string
        end
        % Check if the number of molecules and reactions are the same
        if strcmp(sType,options.type)
            % They match
            set_match.type(end+1) = iSet;
        else
            % They do not match
        end
    end
end

end