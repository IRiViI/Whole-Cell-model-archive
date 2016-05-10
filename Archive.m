% Create archive of whole cell model
%   archive = Archive();
%   archive.add_sets('/home/path/to/sets');

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 10/05/2016

classdef Archive < dynamicprops
    
    properties
        
    end
    
    methods
        
        function this = Archive()
            % Initiate archive
            
            % Properties
            this.addprop('settings');
            this.addprop('set');
            this.addprop('info');
            this.addprop('stoichiometry');
            
            % Path to archive class folder
            tmp = dbstack('-completenames');
            % Find the correct dir (apprently required durint parallel processing...)
            % It seems that there are more directories found when using
            % parfor
            %        (tmpNum)     tmpNum = find(arrayfun(@(x)(~isempty(strfind(tmp(x).file,'Archive.mat'))),1:length(tmp)));
            this.settings.dir = fileparts(tmp.file);
            
            % Add paths
            this.add_paths();
            
            % Set name of archive. This name is used by the function 
            % archive.save_archive.
            this.settings.name = 'new_archive';
            
            % Auto save state. The archive is automatically saved after
            % certain methods when auto_save is true
            this.settings.auto_save = true;
            
        end
        
        function this = add_paths(this)
            % Add paths
            
            addpath(this.settings.dir);
            addpath([this.settings.dir '/methods']);
            addpath([this.settings.dir '/functions']);
            addpath([this.settings.dir '/data_base']);
        end
        
        function this = initiate_WCM(this)
            % Add the libraries of the whole cell model
            
            % Save current directory
            current_directory = cd;
            % Check how many times a new directory is set
            attempt = 1; 
            
            % Try to open while not succeeded
            while true
                
                try % Try
                    % Open data base of the whole cell model
                    [~]=cd(this.settings.WCM);
                    addpath(this.settings.WCM);
                    setWarnings();
                    setPath();
                    
                    % Go back to original directory
                    cd(current_directory);
                    
                    % Break after succes
                    break
                    
                catch % If fails
                    
                    % Check if the directory of the whole cell model is know
                    if ~isfield(this.settings,'WCM')
                        % If no directory is set so far
                        tmp = input('Please give the directry of the whole cell model on local computer:\n','s');
                        this.settings.WCM = tmp;
                    else
                        if attempt > 1 % If it's the second attempt
                            multiAttempt = 'Example: /home/path/to/folder/wholeCell-master\n';
                        else
                            multiAttempt = '';
                        end
                        % If the directory is set but incorrect
                        inputMessage = sprintf(['Whole cell model directory is incorrect.\n' ...
                            'Current directory is \n%s\n'...
                            '%s'...
                            'please enter correct directory:\n'],this.settings.WCM,multiAttempt);
                        tmp = input(inputMessage,'s');
                        this.settings.WCM = tmp;
                    end
                    
                    % Update attempt
                    attempt = attempt + 1;
                    
                end
                
            end
            
            % Update the user
            fprintf('Whole cell model libararies opened\n');
            
        end
        
        function this = add_sets(this, varargin)
            % Add sets and/or simulations to existing sets
            %
            %   archive.add_sets('path/to/simulations/dir')
            %
            
            add_sets_archive(this,varargin{:});
            
            this.save_archive(true);
            
        end
        
        function this = general_processing(this)
            
            %             fprintf('Determine life time:\n');
            %             this.life_time();         % OLD DATED
            %             fprintf('Determine growth rate:\n');
            %             this.growth_rate();         % OLD DATED
            fprintf('Determine trends:\n');
            this.trends();
            
        end
        
        function this = general_information(this)
            % Get some general information about the whole cell model
            
            fprintf('Set labels:\n');
            this.set_labels();
            fprintf('Done\n');
            fprintf('Set stoichiometry:\n');
            this.set_stoichiometry();
            fprintf('Set categories to info list:\n');
            this.extract_from_xls_file('info.Metabolite.ID','info.Metabolite.Category',...
                'mmc4.xlsx','S3G-Metabolites',1,16);
            fprintf('Done\n');
            
        end
        
        function state = load_file(this,varargin)
            % Load file belonging to a simulation. The default state file loaded
            % is "point_compressed". This can be changed by setting
            % 'fileTypeTag' to another type. All fields are loaded unless
            % mentioned otherwise
            %
            % state = archive.load_file(set,simulation)
            %
            %   'fileTypeTag'   - Change the type of state file to load.
            %   'fields'        - Load specific field(s). Mention the field in the form
            %                     of a char array. Example: {'a','b'}.
            
            set = varargin{1};
            simulation = varargin{2};
            
            state = load_state(this,set,simulation,varargin{3:end});
            
        end
        
        function value = extract_data(this, varargin)
            % Extract value of sets or simulations of a specific set
            %
            % value = archive.extract_data({'subfield','subsubfield',...})
            % the fixed field is set, when the target is 'strain' (set).
            %
            %   'set'       - It's possilbe to extract information of
            %                 simulations of a specific set. Set 'set' to
            %                 the set that should be inspected.
            %   'select'    - The coordinates of the value to extract. If the value
            %                 is within in a vector or matrix, specify the x and y
            %                 values Example: list.data(1).field.subfield = rand(3,3)).
            %                 'select',[2,2],... selects the element in the middle.
            %                 Other options are "'select','all'" and "'select','end'".
            
            value = extract_data_archive(this, varargin{:});
            
        end
        
        function this = check_crashed_simulations(this,varargin)
            % Check for crashed simulation
            %
            %   'set'           - Only check one or multiple specific sets.
            %                     Example: 'set',[1,2,10],...
            
            this = crash_check_archive(this,'remove','Off',varargin{:});
            
            this.display_crashed_simulations;
            
        end
        
        function display_crashed_simulations(this)
            % Display reaction that crashed
            
            tSet = length(this.set);
            for iSet = 1:tSet;
                tSim = length(this.set(iSet).simulation);
                for iSim = 1:tSim
                    value=this.set(iSet).simulation(iSim).crash;
                    if value>0;
                        fprintf('crashed: set %d, sim %d\n',iSet,iSim);
                    end
                end
            end
            
        end
        
        function this = life_time(this,varargin)
            % The life time is checked during checking for crashed
            % simulations.
            %
            %   'set'           - Only check one or multiple specific sets.
            %                     Example: 'set',[1,2,10],...
            %   'redo'          - Redo old calculations (Default 'Off')
            
            this = crash_check_archive(this,'remove','Off',varargin{:});
            
        end
        
        function this = remove_crashed_simulations(this,varargin)
            % Remove Crashed simulations
            
            tSet = length(this.set);
            for iSet = 1:tSet;
                tSim = length(this.set(iSet).simulation);
                for iSim = tSim:-1:1 % Backwards in order to keep the order correct after removing sim
                    value=this.set(iSet).simulation(iSim).crash;
                    if value>0;
                        this.set(iSet).simulation(iSim) = [];
                        fprintf('Removed: set %d, sim %d\n',iSet,iSim);
                    end
                end
            end
            fprintf('Completed\n');
            
        end
        
        function this = check_unique_simulations(this)
            % Check uniqueness of seeds within a set
            
            unique_check_archive(this,'remove','Off');
            
            this.display_non_unique_simulations;
        end
        
        function display_non_unique_simulations(this)
            % Display reaction that are not unique
            
            tSet = length(this.set);
            for iSet = 1:tSet;
                tSim = length(this.set(iSet).simulation);
                for iSim = 1:tSim
                    value=this.set(iSet).simulation(iSim).unique;
                    if value==0;
                        fprintf('not unique: set %d, sim %d\n',iSet,iSim);
                    end
                end
            end
            
        end
        
        function this = remove_non_unique_simulations(this)
            % Check and remove non unique seeds within a set
            
            tSet = length(this.set);
            for iSet = 1:tSet;
                tSim = length(this.set(iSet).simulation);
                for iSim = tSim:-1:1
                    value=this.set(iSet).simulation(iSim).unique;
                    if value==0;
                        fprintf('Removed: set %d, sim %d\n',iSet,iSim);
                        this.set(iSet).simulation(iSim) = [];
                    end
                end
            end
            fprintf('Completed\n');
            
        end
        
        
        function plot(this,varargin)
            % Plot simulations
            %
            %   plot(set,{field},compartment)
            %
            %   Examples:
            %   archive.plot(9,{'Metabolite','counts'},50,1)
            %   archive.plot(9,{'mass'},1)
            %
            %   'figureNumber'  - Set the figure number of the plot. Default generates a
            %                     new figure.
            %   'fileTypeTag'   - Target file that should be used for this analysis.
            %                     This function looks for files that contain a certain
            %                     word (or string) in the file name.
            %   'simulation'    - Select simulations to plot (example [1,2,5])
            %   'symbol'        - Set the symbol for the data points (example 'o')
            %   'step'          - Set the step size. Meaning, how many data points
            %                     should be skipped when plotting (example 10)
            %   'glueState'     - Glue togheter the multiple state files when needed.
            %                     'glueState' can either be 'mean' or 'end'. meaning
            %                     the mean or last value of the state file.
            %   'progressState' - Show progress. ('On' or 'Off')
            %   'Color'         - Color of symbol
            
            plot_archive(this,varargin{:})
        end
        
        function cyto_viz(this,varargin)
            % Visualize data
            %
            % archive.cyto_viz(''
            
            % Default settings
            options.file = 'Default_layout.xgmml';
            options.reactID = this.info.MetabolicReaction.ID;
            options.metaboliteID = this.info.Metabolite.ID;
            options.progress = 'On';
            options.Metabolite = 'Off';
            options.MetabolicReaction = 'Off';
            
            % Adjust options
            inputOptions = struct(varargin{:});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % Default values if left flux empty
            tReaction = length(this.info.MetabolicReaction.ID);
            colorData = 0.5*ones(1,tReaction);
            thicknessData = 0.5*ones(1,tReaction);
            
            % Apply flux values
            if strcmp(options.MetabolicReaction,'On')
                if isfield(options,'colorData')
                    colorData = options.colorData;
                end
                if isfield(options,'thicknessData')
                    thicknessData = options.thicknessData;
                end
            end
            
            % Default setting metabolite data
            metaboliteData = [];
            
            % Adjust default setting metabolite data
            if strcmp(options.Metabolite,'On')
                metaboliteData = options.metaboliteData;
            end
            
            add_data_xml(options.file,options.reactID,...
                thicknessData,...
                colorData,...
                'outputFile','temp.xgmml',...
                'progress',options.progress,...
                'metaboliteID',{options.metaboliteID},...
                'metaboliteData',metaboliteData,...
                'metabolicReaction',options.MetabolicReaction);
            
        end
        
        function this = set_labels(this, varargin)
            % Set labels
            %
            %   Optional:
            %
            %   'excelFile'     - Location of excel file
            
            % Default file
            options.excelFile = 'StatePropertyRowColIDs.xlsx';
            
            % Adjust options
            inputOptions = struct(varargin{:});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % Add labels to archive
            this = labels_archive(this, options.excelFile);
            
        end
        
        function this = set_stoichiometry(this, varargin)
            % Get and implement stoichiometry info
            
            % Default file
            options.excelFile = 'mmc4.xlsx';
            options.progress = 'On';
            
            % Adjust options
            inputOptions = struct(varargin{:});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            if strcmp(options.progress,'On')
                fprintf('Read xml file:\n')
            end
            
            % Add stoichiometry info
            this = stoichiometry_archive(this, options.excelFile);
            if strcmp(options.progress,'On')
                fprintf('Done\nProcess stoichiometry information:')
            end
            % Implement stoichiometry structure;
            this = extent_stoichiometry_archive(this, 'progress', options.progress);
            
        end
        
        function compress_state_files(this,varargin)
            % Compress the state files. Construct two files in which the
            % last or mean values of all the parameters are stored for
            % each state file. The files are called 'mean_compressed_state.mat' and
            % 'point_compressed_state.mat'.
            %
            %   'set'           - Select specific set. Set set to "0" to process all
            %                     sets
            %   'simulation'    - Select specific simualtion. Set simulation to "0" to
            %                     process all simulations
            %   'redo'          - redo simultions (true or false). default: false
            
            options.minimal = false;
            
            % Adjust options
            inputOptions = struct(varargin{:});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % Check if whole cell model libraries are opened
            if ~options.minimal
                this.check_WCM_library('ask',true)
            end
            
            % Execute compression
            compress_state_files_archive(this,varargin{:});
            
        end
        
        function [setAverage, simulationvalues] = average_value(this,varargin)
            % Get average value of a set (also gives the values of the
            % simulations)
            %
            % [strainAverage, simulationValues] = average_value(...
            %  set,timePoint,field,subfield,compartment);
            %
            % NOTE: The default uses the 'mean_compressed' file to get the
            % average
            
            [setAverage, simulationvalues] = average_value_archive(this,varargin{:});
            
        end
        
        function std_snap_shot = deviation_snap_shot(this,varargin)
            % It is now quite unuseful because I incorporated it in
            % snap_shot
            
            iSS = varargin{1};
            iSet = varargin{2};
            
            total_snap_shot = [];
            
            tSim = length(this.set(iSet).simulation);
            
            for iSim = 1:tSim
                snap_shot = this.set(iSet).simulation(iSim).snapShot(iSS).values;
                
                total_snap_shot(:,:,:,iSim) = snap_shot;
                
            end
            
            std_snap_shot = nanstd(total_snap_shot,[],4);
            
        end
        
        function this = classification(this, varargin)
            % Classify the simulations and sets
            %
            %   'set'           - Only check one or multiple specific sets.
            %                     Example: 'set',[1,2,10],...
            
            %             try
            this = classification_archive(this, varargin{:});
            %             catch
            %                 warningMessage = sprintf(['Note: this method requires the methods "trends", \n'...
            %                     '"growth_rate" and "life_time" \n'...
            %                     ' to be executed first\n']);
            %                 warning(warningMessage);
            %             end
            
        end
        
        function this = trends(this, varargin)
            % Get the trends of the simulations and sets
            %
            %   'set'           - Only check one or multiple specific sets.
            %                     Example: 'set',[1,2,10],...
            %   'fileTypeTag'   - Select the tag that the folder should contain. The
            %                     default is value is "summary" such that the summary
            %                     files are analyzed.
            %   'fields'        - Select the fields that should be fitted. Default:
            %                     fields = {'proteins','rnas'}. Note: This function
            %                     might not support a certain field. This function
            %                     should be adjusted accordingly if this is the case.
            %   'progressState' - Show progress of this function ('On' or 'Off').
            %                     Default value is 'Off'.
            %   'redo'          - Redo old calculations (Default 'Off')
            
            this = trends_archive(this, varargin{:});
            
            this.save_archive(true)
            
        end
        
        function this = growth_rate(this, varargin)
            % Get the growth rates of the simulations and sets
            %
            %   'set'           - Only check one or multiple specific sets.
            %                     Example: 'set',[1,2,10],...
            %   'progress'      - Display the progress of this function. Default
            %                     setting is 'Off'
            %   'plot'          - Plot the actual data and the fit while running this
            %                     function. Default setting is 'Off'.
            %   'fileTypeTag'   - Target file that should be used for this analysis.
            %                     This function looks for files that contain a certain
            %                     word (or string) in the file name.
            %   'redo'          - Redo old calculations (Default 'Off')
            
            this = growth_rate_archive(this, varargin{:});
            
        end
        
        function this = snap_shot(this, varargin)
            % Make a snap shot of every simulation at a certain time point
            % New snap shots are stored additionally to the pre-existing
            % snap shots. The settings of the snap shots can be found in
            % the achive.settings.snap_shot section.
            %
            %   'compartment'   - Set the compartment number. Default = 1
            %   'fileTypeTag'   - Make a snap shot of a certain file. Default:
            %                     'mean_compressed'
            %
            % Example: archive.snap_shot(2500,'Metabolite','counts')
            
            this = snap_shot_archive(this, varargin{:});
            
            this.save_archive(true);
            
        end
        
        function this = remove_snap_shot(this, varargin)
            % Remove a or all snap shots of archive.
            %
            % archive.remove_snap_shot(0) removes all snap shots
            % archive.remove_snap_shot(1:5) removes snap shots 1:5
            
            if nargin > 1
                target = varargin{1};
            else
                error('Not enough input arguments');
            end
            
            if target == 0
                this.settings.snapShot(:) = [];
            else
                this.settings.snapShot(target) = [];
            end
            
            % Number of sets
            tSet = length(this.set);
            
            for iSet = 1:tSet
                
                % Number of simulations
                tSimulation = length(this.set(iSet).simulation);
                
                % Remove snap shots from simulations
                for iSimulation = 1:tSimulation
                    if target == 0
                        this.set(iSet).simulation(iSimulation).snapShot = [];
                    else
                        this.set(iSet).simulation(iSimulation).snapShot(target) = [];
                    end
                end
                
                % Remove snap shot from set
                try
                    if target == 0
                        this.set(iSet).snapShot = [];
                    else
                        this.set(iSet).snapShot(target) = [];
                    end
                catch
                end
                
            end
            
        end
        
        function varargout = dendrogram(this, varargin)
            % Method archive.snap_shot required for this method.
            %
            % Make a dendogram of snap shot
            %
            % [z,selection] = archive.dendrogram(snap_shot_number,...)
            % Where "z" is the z matrix created for the construction of the
            % dendrogram and "selection" includes the reactions or
            % molecules which are used in construction of the dendrogram.
            %
            % There are many options (see dendrogram_archive.m). Including:
            %
            %   'compare'       - Make a dendrogram using the average values of the
            %                     set ('average') or of individual simulations
            %                     ('individual').
            %   'reference'     - Scale the values are skilled using a set
            %                     ('average') or a simulation ('individual'). Default:
            %                     'average'
            %  'referenceNumber'- The number of the set or simulation which should
            %                     be used as reference. (yeah, at the moment selecting
            %                     a simulation as reference is not straight forward if
            %                     it's not within the first strain). Default: 1
            %   'selection'     - Tresshold value: "values" obtained after the scalling
            %                     should be greater than this value.
            %   'include'       - Reactions or molecules to include (default 'all')
            %   'set'           - Sets to include (default 'All')
            %   'colortreshold' - The standard "colortreshold" option of the function
            %                     dendrogram.m
            %
            %   Examples:
            %   z = archive.dendrogram(1,'include','all',...
            %       'compare','average','set',[1,2,3,4,5,10],'selection',1);
            %   z = archive.dendrogram(1,'set',archive.extract_data({'class'})==1);
            
            % Snap shot number
            snap_shot_number = varargin{1};
            
            % Make dendrogram
            [z,selection,D]=dendrogram_archive(this, snap_shot_number, varargin{2:end});
            
            % varargout could probably be directly linked but a well.
            varargout{1} = z;
            varargout{2} = selection;
            varargout{3} = D;
        end
        
        function [this,track_clusters,cluster_list] = clustering(this, varargin)
            % method archive.dendrogram required for this method.
            %
            % [track_clusters,cluster_list,archive] = find_clusters_archive(z,...
            %   nClusters,sClusters,set,selection,snap_shot_number)
            %
            % Cluster the sets according to the z matrix generated with the
            % archive.dendrogram method. This function requires: the z
            % matrix, the number of clusters that you think are present in
            % the dendrogram, the minimal size of the clusters, the sets
            % which are involved in the creation of the dendrogram, the
            % molecules/reactions whicha are used for the construction of
            % the z matrix (given as "selection" by the archive.dendrogram
            % method), and the number of the snap shot used for creating
            % the dendrogram
            %
            % Examples:
            % [z,selection] = archive.dendrogram(1)
            % [track_clusters,cluster_list,archive] = archive.clustering(z,...
            % 3,20,[1:length(archive.set)],selection,1)
            %
            % [z,selection] = archive.dendrogram(1,'set',[1,1,1,0,1])
            % [track_clusters,cluster_list,archive] = archive.clustering(z,...
            % 3,20,[1,1,1,0,1],selection,1)
            
            z = varargin{1};
            nClusters = varargin{2};        % How many clusters you want
            sClusters = varargin{3};        % The minimal size of the clusters
            reactions_or_molecules_used = varargin{4}; % The reactions or molecules involved in the construction of z ('selection' of archive.dendrogram)
            snap_shot_number = varargin{5}; % Number of the snap shot used in order to construct the z matrix
            
            %             if
            %             set = varargin{4};              % The sets involved in the construction of the z matrix
            %             end
            
            [track_clusters,cluster_list,this] = find_clusters_archive(z,...
                nClusters,sClusters,this,reactions_or_molecules_used,...
                snap_shot_number);
        end
        
        function histogram(this,varargin)
            % archive.histogram(snap_shot_number, element_coordinates_matrix)
            
            
            options.snap_shot = varargin{1};
            options.element = varargin{2};
            
            options.fields = fields(this.settings.classification);
            
            % Number of classes
            tClass = length(options.fields);
            
            % Initiate individual histogram structure
            for iClass = 1:tClass
                individual(iClass).value = [];
            end
            
            % Initiate total histogram structure
            totalValue = [];
            
            % Total number of sets
            tSet = length(this.set);
            
            for iSet = 1:tSet
                
                % Total number of simulations
                tSim = length(this.set(iSet).simulation);
                
                for iSim = 1:tSim
                    
                    % Get class
                    class = this.set(iSet).simulation(iSim).class;
                    
                    % Get value
                    value = this.set(iSet).simulation(iSim).snapShot(options.snap_shot).values(options.element);
                    
                    % Save values
                    totalValue(end+1) = value;
                    individual(class+1).value(end+1) = value;
                    
                end
                
            end
            
            % Get settings
            figure()
            temp = histogram(totalValue);
            binWidth = get(temp,'binWidth');
            clf;
            
            % Add histograms
            hold on
            for iClass = 1:tClass
                histogram(individual(iClass).value,'binWidth',binWidth);
            end
            hold off
            
        end
        
        function this = extract_from_xls_file(this,varargin)
            %Extract data of a xls file and add it to the archive.
            %
            % archive = xls_to_archive(archive,'field1.field2.linker',...
            %   'field1.field2.data','path/to/xls/file','xls_sheet_name',linker_column,
            %   data_column);
            %
            %  In this file are the values in the column "data_column" of the
            %  "xls_sheet_name" in the "path/to/xls/file" added to the "archive" in the
            %  subfield "archive.field1.field2.data". The values are inserted in the order
            %  according to the linkers.
            %
            % Example: element(1,42) is 'ATP' in the data sheet and
            % archive.field1.field2.linker(33) or archive.field1.field2(33).linker is also
            % 'ATP' then: archive.field1.field2.data(33) and archive.field1.field2(33).data
            % are the value of element(4,42) respectively in the case where data_column
            % is 4.
            %
            % Note:
            % There are two different ways that the archive is structured. In one case is
            % are all the values saved as one array under the subfield
            % (archive.field.subfield = array). In the other case, are all the values
            % saved under different entries of the field (archive.field(1).subfield =
            % array(1), archive.field(2).subfield = array(2),...). Both cases are taken
            % into account
            
            this = xls_to_archive(this,varargin{:});
            
        end
        
        function this = insert_set(this,varargin)
            % Insert an existing set to another location
            % archive.instert_set(old_location,new_location);
            
            % The target and insert
            target = varargin{1};
            insert = varargin{2};
            
            % Number of sets
            tSet = length(this.set);
            
            % Make new order list
            order = 1:tSet;
            order(insert:end) = order(insert:end) -1;
            if target > insert
                if target~=tSet
                    order(target+1:end) = order(target+1:end) + 1;
                end
            elseif target < insert
                if target~=tSet
                    order(target:end) = order(target:end) + 1;
                end
            end
            order(insert) = target;
            
            % Just checking uniqueness (Quite unnecessary, but a well)
            check_order = unique(order);
            if length(check_order)~=length(order)
                error('Something went wrong while ordering sets!');
            end
            
            % Apply new order
            this.set = this.set(order);
            
        end
        
        function health_check(this, varargin)
            % Check the health of file by trying to load it.
            % archive.health_check(target)
            % target should be the name of the file.
            %
            % Example:
            %   archive.health_check('mean_compressed_state.mat')
            %
            %   'compress'      - Compress files if unhealhty file is
            %                     found. Default 'Off'
            
            % Mandatory inputs
            target = varargin{1};
            
            % Default settings
            options.compress = 'Off';
            
            % Adjust options
            inputOptions = struct(varargin{2:end});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % Progress spacer
            fprintf('Progress:\n      ');
            
            % Number of sets
            tSet = length(this.set);
            
            % For every set
            for iSet = 1:tSet
                
                % Progress
                display_progress(iSet,tSet);
                
                % Total number of simulation
                tSim = length(this.set(iSet).simulation);
                
                % For every simulation
                for iSim = 1:tSim
                    
                    % Try to load file
                    try
                        this.load_file(iSet,iSim,'fileTypeTag',target);
                    catch
                        % If fails:
                        fprintf('Unhealthy: Set %d, Simulation %d\n',iSet,iSim);
                        
                        % Recompress files if indicated
                        if strcmp(options.compress,'On')
                            fprintf('Compressing:\n')
                            this.compress_state_files('set',iSet,'simulation',iSim,'redo','On');
                            try
                                % Try loading file again
                                this.load_file(iSet,iSim,'fileTypeTag',target);
                                fprintf('Works correctly now!\n')
                            catch
                                fprintf('File does not function\n')
                            end
                        end
                        
                        fprintf('      ');
                        
                    end
                    
                end
                
            end
            
        end
        
        function save_figure(this,varargin)
            % Save current figure
            % archive.save_figure('file_name')
            
            save_figure_archive(this,varargin{:})
        end
        
        function save_archive(this,varargin)
            % Save archive
            % Change the value archive.settings.name in order to change the
            % auto save name.
            %
            %   Examples:
            %   archive.save_archive
            %       Always save simulation
            %
            %   archive.save_archive(true)
            %       Only save if auto save is true
            %
            
            % Set auto save
            if nargin > 1
                auto = varargin{1};
            else
                auto = false;
            end
            
            % Save archive if desired
            if this.settings.auto_save || (~auto)
                
                % Change structure name archive
                archive = this;
                
                % Archive name
                name = this.settings.name;
                
                % Check for additions
                if auto
                    variable = 'Auto_save_';
                else
                    variable = '';
                end
                
                % Update name
                full_name = [variable name '.mat'];
                
                % Save archive
                save(full_name,'archive');
                
                % Inform user
                fprintf('Saved: %s\n', full_name);
            end
            
        end
        
        function sets = sets(this)
            % tSets = archive.sets
            %
            % Number of sets in archive
            
            sets = length(this.set);
            
        end
        
        function sims = simulations(this, varargin)
            % tSimulation = archive.simulations(set)
            %
            % Number of sims in sets "set"
            % if set = 0, then all simulations are counted in all sets.
            
            % Process input arguments
            if nargin == 1
                set = 0;
            elseif nargin == 2
                set = varargin{1};
            else
                warning('Too many input arguments');
            end
            
            % Total number of sets
            tSet = this.sets;
            
            if set == 0
                lSet = 1:tSet;
            else
                lSet = set;
            end
            
            % Initiate sims counter
            sims = 0;
                
            % For every set
            for iSet = lSet
                
                % Add sims
                sims = sims + length(this.set(iSet).simulation);
                
            end
                
        end
        
        function this = add_info(this, varargin)
            % Add a new info structure to archive
            %
            % archive.add_info('path/to/infoFile')
            %   or
            % archive.add_info(set) % where set is a numeric
            %
            %   'note'      - Add note to info file
            
            % File dir of info.mat file
            if isstr(varargin{1})
                dir = varargin{1};
                load(dir); % Load info structure
            else isnumeric(varargin{1})
                % Look for info.mat files in folders
                tSim = length(this.set(varargin{1}).simulation);
                for iSim = 1:tSim
                    folder = this.set(varargin{1}).simulation(1).folder;
                    dir = [folder '/' 'info.mat'];
                    try
                        load(dir); % Load info structure
                        break % Stop if info.mat structure found
                    catch
                        % If it's the last simulation folder
                        if iSim == tSim
                            error('No info.mat found in simulation folders');
                        end
                    end
                end
            end
            
            
            % Default settings
            options.note = 'new info structure';
            
            % Check if the basic info structure is present
            if isempty(this.info)
                error('Please include the standard information default first. archive.general_information')
            end
            
            % Optional parameters
            % Adjust options
            inputOptions = struct(varargin{2:end});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % The fields and number of fields of preexisting info structure
            prefInfo = fields(this.info(1));
            pretInfo = length(prefInfo);
            
            % Add preexisting fields to new info structures
            info = struct();
            info.note = options.note;
            for preiInfo = 1:pretInfo
                % Check if field exist
                if ~isfield(info,prefInfo{preiInfo})
                    % Add field if it does not exist
                    info.(prefInfo{preiInfo}) = struct;
                end
            end
            
            % The fields and number of fields of preexisting info structure
            fInfo = fields(info);
            tInfo = length(fInfo);
            
            % Add new fields to preexisting info structure
            for iInfo = 1:tInfo
                % Check if field exists
                if ~isfield(this.info,fInfo{iInfo})
                    % Add field if it does not exist
                    this.info.(fInfo{iInfo}) = struct;
                end
            end
            
            % The fields and number of newfields info structures
            newfInfo = fields(this.info(1));
            newtInfo = length(newfInfo);
            
            % Order the new info structure into the same order as
            % preexisting info structures
            for newiInfo = 1:newtInfo
                shuffledinfo.(newfInfo{newiInfo}) = info.(newfInfo{newiInfo});
            end
            
            % Add info structure
            this.info(end+1) = shuffledinfo;
            
        end
        
        function out = get_files(this, varargin)
            % archive.get_files(fileTypeTag,set)
            %
            % 'simulation'  - Specify specific simulation
            %   'warning'   -'On' or 'Off', supress warning statements of this function.
            %                 Default: 'On'
            %   'multi'     -'On' or 'Off', multiple files are allowed in a folder and
            %                 will be added to the archive. Default: 'On'
            
            out = get_files_archive(this, varargin{:});
            
        end
        
        function this = inspect_state_files(this,varargin)
            % Obtain the time range of every individual "state-x" file
            %
            %   'set'       - Set which set to use
            %   'simulation'- Set which set to use
            
            % Default settings
            options.simulation = 1;
            options.set = 1;
            
            % Adjust options
            inputOptions = struct(varargin{:});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % Get path to the state files
            out = this.get_files('state-',options.set,'simulation',options.simulation);
            
            % Total number of state files
            tState = length(out);
            
            % Initate state file information structure
            info = struct();
            
            % Progress
            fprintf('Progress:\n     ');
            
            % For every state file
            for iState = 1:tState
                
                % Progress
                display_progress(iState,tState)
                
                % Split string to get the name of the state file
                tmp = strsplit_archive(out{iState},'/'); % Remove "/"
                tmp = strsplit_archive(tmp{end},'.'); % Remove "."
                tmp = strrep(tmp, '-', '_'); % Replace "-"
                tmp = strrep(tmp, ' ', '_'); % Replace " "
                niState = tmp{1}; % The name of the state file
                
                % Load the state file
                state = load(out{iState},'Time');
                
                % Get values
                info.(niState).Time.min = min(squeeze(state.Time.values));
                info.(niState).Time.max = max(squeeze(state.Time.values));
            
            end
            
            this.info.state = info;
        end
        
        function find_time(this,varargin)
            % Look which state file types might include the desired time
            % point
            
            % Settings
            options.time = varargin{1}; % time point
            options.info = 1; % The state info number
            
            % Total number of state files
            fState = fields(this.info.state);
            tState = length(fState);
            
            % For every type of state file
            for iState = 1:tState
                
                % Get the time info structure
                info = this.info(options.info).state.(fState{iState}).Time; 
                
                % Check if the time point is between the two extreems
                if info.min <= options.time && options.time <= info.max
                    
                    fprintf('The time point %d is found in state file type %s\n',...
                        options.time,fState{iState});
                    
                end
                
            end
            
        end
        
        function out = count_simulation_folder_content(this)
            % Counter how many simulations are found in the folders
            
            % Number of sets
            tSet = this.sets;
            
            % Inititate out structure
            out = struct();
            
            % For every set
            for iSet = 1:tSet
                
                % Number of simulations in set
                tSim = this.simulations(iSet);
                
                % For every simulation of set
                for iSim = 1:tSim
                    
                    % Get the folder path
                    folder = this.set(iSet).simulation(iSim).folder;
                    
                    % Cut folder path
                    element = strsplit_archive(folder,'/');
                    
                    % total number of elements
                    tElement = length(element);
                    
                    % For every element length
                    for iElement = 1:tElement
                        
                        % Only process if the element is not empty
                        if ~strcmp(element(iElement),'')
                            
                            % Create path
                            path = [element(1:iElement)];
                            path = path(~strcmp(path,''));
                            path = strrep(path,' ','_');
                            tPath = length(path);
                            for iPath = 1:tPath
                               path{iPath} =  ['folder_' path{iPath}];
                            end
                            
                            % Path to counter
                            path_counter = {path{:},'simulations'};
                            
                            % Update counter
                            try
                                tmp = getfield(out,path_counter{:});
                                out = setfield(out,path_counter{:},tmp+1);
                            catch
                                out = setfield(out,path_counter{:},1);
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function out = infocmp(this,varargin)
            % strcmp like function specializid for archive which uses
            % strfind instead of strcmp. 
            %
            %   out = infocmp(target,{'field1','field2',..})
            %   
            %   Example:
            %       out = infocmp('DNA',{'info','ProteinComplex','name'});
            %       out = 
            %           [1 1 1 ..]
            
            % Get structure
            if nargin > 3
                field = get_field_archive(this, {varargin{2}{1:end-1}},varargin{3});
            else
                field = get_field_archive(this, {varargin{2}{1:end-1}});
            end
            tmp = {field.(varargin{2}{end})};
            
            % Compare the tag in the structure
            out = arrayfun(@(x)(~isempty(strfind(tmp{x},varargin{1}))),1:length(tmp));
            
        end
        
        function cell_cycle(this,varargin)
            % Visualize the progress of the cell
            %
            %   'set'       - set number
            %  'simulation' - simulation number
            
            % Check if whole cell model libraries are loaded; load
            % libraries if missing
            this.check_WCM_library('initiate',true);
            
            % Execute analysis
            cell_cycle_archive(this, varargin{:});
            
        end
        
        function check_WCM_library(this, varargin)
            % Check if the library is added
            
            % Default settings
            options.initiate = false;
            options.ask = false;
            
            % Adjust options
            inputOptions = struct(varargin{1:end});
            fieldNames = fieldnames(inputOptions);
            for i = 1:size(fieldNames,1)
                options.(fieldNames{i}) = inputOptions.(fieldNames{i});
            end
            
            % Get folders
            tmp = path;
            
            % Find WCM library part
            tmp = strfind(path,'/lib/absolutepath');
            
            % Process state
            if isempty(tmp)
                
                % Ask if they want to load the libraries
                if options.ask
                    while true
                        answer = input(['Libraries of the whole cell model are missing\n'...
                            'Would you like to load whole cell model libraries?\n'...
                            '(yes/no)\n'],'s');
                        if strcmpi(answer,'yes')
                            options.initiate = true;
                            break
                        elseif strcmpi(answer,'no')
                            options.initiate = false;
                            break
                        else
                            fprintf(' -- INVALID INPUT -- \n')
                        end
                    end
                end
                
                % Load libraries
                if options.initiate
                    this.initiate_WCM;
                else
                    warningMessage = sprintf('It''s adviced to inititate whole cell model library first.\n Execute: archive.initiate_WCM.');
                    warning(warningMessage);
                end
            end
            
        end
        
        function this = template(this,varargin)
            % This is template example for adding new function to the
            % archive structure.
            % archive.('setting1',setting1,'setting2',setting2,...)
            
            % Optionals:
%             this.check_WCM_library;
            
            % Execute function
            this = template_archive(this,varargin{:});
            
            % Optionals:
%             this.save_archive(true);
            
        end
    end
    
end
