

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function cmean_state_file_archive(archive, varargin)
%% Settings

% Default settings
options.set = 0;        	% Default: analyze all sets
options.simulation = 0; 	% Default: analyze all simulations
options.progress = true;    % Default: show progress
options.fileTypeTag = 'mean_compressed_state.mat';
options.dirTemp = [archive.settings.dir '/' 'temp'];
options.dirOut = [archive.settings.dir '/' 'output'];
options.uuid = char(java.util.UUID.randomUUID); % Random code for temp files

% Change settings
tmpOptions = struct(varargin{1:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

%% Process settings

% Determine which sets of the archive should be analyzed
if options.set == 0
    tSet = archive.sets;    % Number of sets to analyze
    lSet = 1:tSet;          % List of sets to analyze
else
    lSet = options.set;     % List of sets to analyze
    tSet = length(lSet);    % Number of sets to analyze
end

% Determine number of simulations to analyze (required for progess report)
if options.simulation == 0
    options.ttSimulation = archive.simulations(lSet);
else
    options.ttSimulation = tSet * length(options.simulation);
end

% Output name
if ~isfield(options,'output')
    tmp = strsplit(options.fileTypeTag,'.');
    options.output = tmp{1};
end

% Create folders if it doesn't exist yet
[~,~,~] = mkdir(options.dirTemp);
[~,~,~] = mkdir(options.dirOut);

% Load file to get field types
state = archive.load_file(lSet(1),1,'fileTypeTag',options.fileTypeTag);

%% Run code

% Initiate simulation counter
options.iiSimulation = 0; % Simulation counter
options.progressCount = 0;

% Display progress
if options.progress
    fprintf('Progress:\n     ');
    display_progress(options.iiSimulation, options.ttSimulation)
end
    
% Make mean and std state files of the following sets
for iSet = lSet
    
    clear mean_state std_state
    
    % Make list of sets to analyze
    if options.simulation == 0
        tSimulation = archive.simulations(iSet);    % Number of simulations to analyze
        lSimulation = 1:tSimulation;                % List of simulations to analyze
    else
        lSimulation = options.simulation;                  % List of simulations to analyze
        tSimulation = length(lSimulation);          % Number of simulations to analyze
    end
    
    % Name temporary total_state folder
    options.path_tmp_total_state = [options.dirTemp '/' 'mean_std_state' '_' num2str(iSet) '_' options.uuid];
    [~,~,~] = mkdir(options.path_tmp_total_state);
    
    % Get field information
    lField = fields(state);
    tField = length(lField);
    
    for iField = 1:tField
        
        % current field
        field = lField{iField};
        fprintf([field '      '])
        lSubField = fields(state.(field));
        tSubField = length(lSubField);
        
        for iSubField = 1:tSubField
            
            location = 0;
            
            clear mean_out std_out
            
            while true
                
                location = location + 1;
                
                clear total_state
                
                % Current subfield(s)
                subfield = lSubField{iSubField};
                
                % Create total state file
                try
                    total_state = create_total_state(archive,iSet,lSimulation,field,subfield, location, options);
                catch
                    break
                end
                
                % Resize total_state matrix
                total_state = resize(archive, total_state, field, subfield, options);
                
                % Creat total mean and std state file
                [mean_out(:,location,:,:), std_out(:,location,:,:)]  = create_mean_std_state(archive,total_state,options);

                % Add subfield to mean and std state structures
                try
                    mean_state.(field).(subfield) = mean_out.(field).(subfield);
                    std_state.(field).(subfield) = std_out.(field).(subfield);
                catch
                    % The dimensions of Polypeptide change and that causes
                    % problems
                    continue
                end
                
            end
            
            
            % Update progress counter
            options.progressCount = options.progressCount + 1/tSubField/tField/tSet;
            
            % Progress
            display_progress(options.progressCount, 1)
            
        end
        
    % Save structures
    save([options.dirOut '/' 'mean_' options.output '_' num2str(iSet)], '-struct', 'mean_state','-append');
    save([options.dirOut '/' 'std_' options.output '_' num2str(iSet)], '-struct', 'std_state','-append');
    
    clear mean_state std_state
    
    end
    
    % Remove temporary folder
%     rmdir(options.path_tmp_total_state,'s');
    
end

% Display progress
if options.progress
    fprintf('Done\n');
end

end

function total_state = create_total_state(archive,iSet,lSimulation, field, subfield, location, options,varargin)
% Create a total_state structure and save it as a temporary file in temp
% folder

% total number of simulations to handle
tSimulation = length(lSimulation);

% Sim counter
iSim = 0;

for iSimulation = lSimulation
    
    % Folder state file
    folder = archive.set(iSet).simulation(iSimulation).folder;
    
    % Path to state file
    pathFile = [folder '/' options.fileTypeTag];
    
    % Load file
    state = load(pathFile,field);
    
    % For to large structures
    iSim = iSim + 1;
    total_state.(field)(iSim).(subfield) = state.(field).(subfield)(:,location,:,:);
    
end

end

function total_state = resize(archive, total_state, field, subfield, options)

tSimulation = length(total_state.(field));
for iSimulation = 1:tSimulation
    lSize(iSimulation,:) = size(total_state.(field)(iSimulation).(subfield));
end
minSize = min(lSize);


fast = true;
if max(strcmp(field,{'Metabolite','ProteinComplex','ProteinMonomer','Rna'}))
   fast = false; 
end

if fast
    
    for iSimulation = 1:tSimulation
        
        total_state.(field)(iSimulation).(subfield) = ...
            total_state.(field)(iSimulation).(subfield)(ones(minSize));
        
    end
    
else
    
    [~,~,~] = mkdir([options.dirTemp '/' options.uuid]);
    
    for iSimulation = 1:tSimulation
        
        fileDir = [options.dirTemp '/' options.uuid '/' 'matrix' num2str(iSimulation)];
        
        matrix = total_state.(field)(iSimulation).(subfield);
        
        save(fileDir,'matrix');
        
    end
    
    clear total_state
    
    for iSimulation = 1:tSimulation
        
        fileDir = [options.dirTemp '/' options.uuid '/' 'matrix' num2str(iSimulation)];
        
        load(fileDir,'matrix');
        
        matrix = matrix(ones(minSize));
        
        save(fileDir,'matrix');
        
    end
    
    for iSimulation = 1:tSimulation
        
        fileDir = [options.dirTemp '/' options.uuid '/' 'matrix' num2str(iSimulation)];
        
        load(fileDir,'matrix');
        
        total_state.(field)(iSimulation).(subfield) = matrix;
        
    end
    
end


end

function add_state_to_total_state(archive,iSet,iSimulation,options)

% Folder state file
folder = archive.set(iSet).simulation(iSimulation).folder;

% Path to state file
pathFile = [folder '/' options.fileTypeTag];

% Load file
state = load(pathFile);

% Get field information
lField = fields(state);
tField = length(lField);

for iField = 1:tField
    
    % current field
    field = lField{iField};
    
    lSubField = fields(state.(field));
    tSubField = length(lSubField);
    
    for iSubField = 1:tSubField
        
        clear total_state
        
        % Current subfield(s)
        subfield = lSubField{iSubField};
        
        % File
        file = [options.path_tmp_total_state '/' field];
        
        try
        tmp = whos('-file',file);
        catch
            clear tmp 
        tmp.name = '';
        end
    
        % Update or create total_state file of this field and subfield
        if sum(strcmp({tmp.name},subfield))>0
            total_state.(field) = load(file, subfield);
            total_state = add_state(total_state,state);
        else
            total_state.(field).(subfield) = state.(field).(subfield);
        end
        
        % Determine save structure
        saveStructure = total_state.(field);
        
        % Save structure
        try
            save(file, '-struct', 'saveStructure','-append');
        catch
            save(file, '-struct', 'saveStructure');
        end
        
        
    end
    
end


end

function total_state = add_state(total_state,state)

% Get all fields of total_state
lField = fields(total_state);
tField = length(lField);

% For every field
for iField = 1:tField
    
    % Get field name of field
    field = lField{iField};
    
    % Get all sub fields of total state
    lSubField = fields(total_state.(field));
    tSubField = length(lSubField);
    
    % For every subfield
    for iSubField = 1:tSubField
        
        % Get subfield name
        subfield = lSubField{iSubField};
        
        % Process if field exist in total_state (some fields are not
        % processed because of dimension dismatch)
        if isfield(total_state.(field),(subfield))
            
            if ~isempty(total_state.(field).(subfield))
                
                % Look at the size of the state files
                newElement = size(state.(field).(subfield));
                totalElement = size(total_state.(field).(subfield));
                
                % Get dimension of new state
                newDimension = length(newElement);
                
                % See which structure has the shortest time
                minTime = min(newElement(newDimension),totalElement(newDimension));
                
                % Get elements
                element = newElement;
                element(newDimension) = minTime;
                
                try
                    % Add to total_state
                    if newDimension == 3
                        total_state.(field).(subfield) = cat(newDimension + 1,...
                            total_state.(field).(subfield)(:,:,1:minTime,:),...
                            state.(field).(subfield)(:,:,1:minTime));
                    elseif newDimension == 2
                        total_state.(field).(subfield) = cat(newDimension + 1,...
                            total_state.(field).(subfield)(:,1:minTime,:),...
                            state.(field).(subfield)(:,1:minTime));
                    elseif newDimension == 1
                        total_state.(field).(subfield) = cat(newDimension + 1,...
                            total_state.(field).(subfield)(1:minTime,:),...
                            state.(field).(subfield)(1:minTime));
                    else
                        error('dimension is different')
                    end
                catch
                    % Remove fields if it doesn't work (easy and dirty way)
                    total_state.(field) = rmfield(total_state.(field),(subfield));
                end
                
            end
            
        end
        
    end
    
end
end

function [mean_state, std_state] = create_mean_std_state(archive,total_state, options)

% Default output
mean_state = struct;
std_state = struct;

% Get field information
lField = fields(total_state);
tField = length(lField);

for iField = 1:tField
    
    % Current field
    field = lField{iField};
    
    lSubField = fields(total_state.(field));
    tSubField = length(lSubField);
    
    for iSubField = 1:tSubField
        
        % Current subfield
        subfield = lSubField{iSubField};

        % Determine mean and std states of field
        [mean_out, std_out]  = determine_mean_std(total_state);
        
        % Add subfield to mean and std state structures
        mean_state.(field).(subfield) = mean_out.(field).(subfield);
        std_state.(field).(subfield) = std_out.(field).(subfield);
        
    end

end
end

function [mean_state, std_state]  = determine_mean_std(total_state)

% Get all fields
lField = fields(total_state);
tField = length(lField);

% For every field
for iField = 1:tField
    
    % Get field name of field
    field = lField{iField};
    
    % Get all sub fields
    lSubField = fields(total_state.(field));
    tSubField = length(lSubField);
    
    % For every subfield
    for iSubField = 1:tSubField
        
        subfield = lSubField{iSubField};
        
        % Get dimension of structure
        dimension = length(size(total_state.(field)(1).(subfield)));
        
%         matrix = cat(dimension+1,total_state.(field).(subfield));
%         
%         mean_state.(field).(subfield) = mean(matrix ,dimension+1);
%         std_state.(field).(subfield) = std(matrix,[],dimension+1);
         
        mean_state.(field).(subfield) = mean(cat(dimension+1,total_state.(field).(subfield)) ,dimension+1);
        std_state.(field).(subfield) = std(cat(dimension+1,total_state.(field).(subfield)),[],dimension+1);
        
%         if max(strcmp(field,{'MetabolicReaction','Metabolite','ProteinComplex','ProteinMonomer','Rna','Polypepeptide'})) == 1 && tSimulation > 5
%             
%             % Get dimension of structure
%             dimension = length(size(total_state.(field).(subfield)));
%             
%             % mean and std targeting the last dimension
%             mean_state.(field).(subfield) = mean(total_state.(field).(subfield),dimension(end));
%             std_state.(field).(subfield) = std(total_state.(field).(subfield),[],dimension(end));
%             
%         elseif strcmp(type,'large')
%             
%             % Get dimension of structure
%             dimension = length(size(total_state.(field)(1).(subfield)));
%             
%             mean_state.(field).(subfield) = mean(cat(dimension+1,total_state.(field).(subfield)),dimension+1);
%             std_state.(field).(subfield) = std(cat(dimension+1,total_state.(field).(subfield)),[],dimension+1);
%         else
%             error('unknown type');
%         end
        
    end
    
end

end

function [mean_state, std_state] = process_state_files(archive,lSet,lSimulation,options)
% Combine all the simulations of the simulations and sets in lSet and
% lSimulations list within the archive object
% Load state in order to get field information
state = archive.load_file(lSet(1),lSimulation(1),...
    'fileTypeTag',options.fileTypeTag);

% Get field information
lField = fields(state);
tField = length(lField);

% Counter
options.iiField = 0;

for iField = 1:tField;
    
    % current field
    field = lField{iField};
    
    lSubField = fields(state.(field));
    tSubField = length(lSubField);
    
    for iSubField = 1:tSubField
        
        % Current subfield(s)
        subfield = lSubField{iSubField};
        
        % Get total_state
        total_state = combine_state_files(lSet,lSimulation, archive, field, subfield, options);
        
        % Determine the mean and std of the total state file
        [mean_out, std_out]  = determine_mean_std(total_state);
        
        % Update mean and std state
        mean_state.(field).(subfield) = mean_out.(field).(subfield);
        std_state.(field).(subfield) = std_out.(field).(subfield);
        
        % Update simulation counter
        options.iiField = options.iiField + 1/tSubField;
        
        % Display progress
        if options.progress
            display_progress(options.iiField, tField)
        end
        
    end
    
end

end


function total_state = combine_state_files(lSet,lSimulation, archive, field, subfield, options)

% Create one structure of all state files
for iSet = lSet
    
    for iSimulation = lSimulation
        
        % Folder
        folder = archive.set(iSet).simulation(iSimulation).folder;
        
        % Path to state file
        pathFile = [folder '/' options.fileTypeTag];
        
        % Load file
        tmpState = load(pathFile,field);
        
        % Select subfield
        state.(field).(subfield) = tmpState.(field).(subfield);
        
        % Clear tmp State file (in order to save memory)
        clear tmpState;
        
        % Add state file to structure
        if exist('total_state','var')
            total_state = add_state(total_state, state);
        else
            total_state = state;
        end
        
    end
    
end

end