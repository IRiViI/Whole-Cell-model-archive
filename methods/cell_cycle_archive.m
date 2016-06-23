% Example
%   cell_cycle_archive(archive,'set',1,'simulation',1)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function varargout = cell_cycle_archive(varargin)
%% ----------------- Preparation processes ----------------- %

% Get archive
archive = varargin{1};

% Default settings
options.fileTypeTag = 'point_compressed_state.mat'; % State file name
options.alpha = 0.25;                               % Opaqueness of regions
options.new = true;                                 % New figure
options.info = 1;                                   % The information structure that should be used in this analysis
options.simulation = 1;                             % Default simulation
options.DnaA.name = {'DnaA-ADP',... % Investigate these enzymes on the chromosome
    'DnaA-ATP',...
    'DnaA-ATP 2mer',...
    'DnaA-ATP 3mer',...
    'DnaA-ATP 4mer',...
    'DnaA-ATP 5mer',...
    'DnaA-ATP 6mer',...
    'DnaA-ATP 7mer',...
    };
options.DnaA.useName = {'DnaA-ATP 2mer',... % The proteins that participate in DnaA evolution
    'DnaA-ATP 3mer',...
    'DnaA-ATP 4mer',...
    'DnaA-ATP 5mer',...
    'DnaA-ATP 6mer',...
    'DnaA-ATP 7mer',...
    };
options.segregation.input = {'ProteinMonomer','CobQ/CobB/MinD/ParA nucleotide binding domain';...
    'ProteinComplex','mraZ protein';...
    'ProteinMonomer','GTP-binding protein Era';...
    'ProteinMonomer','GTPase1 Obg'};

% Process inputs
% mandatory: 'set'
tmpOptions = struct(varargin{2:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

if ~isfield(options,'set')
    error('Please specify the set to analyze (''set'',set)');
end

% Folder
folder = archive.set(options.set).simulation(options.simulation).folder;

%% ----------------- Processes information ----------------- %

% Find DnaA protein numbers in info structure
protein = options.DnaA.name;
tProtein = length(protein);
for iProtein = 1:tProtein
    % Find the DnaA enzymes
    iiProtein = find(archive.infocmp(protein{iProtein},{'info','ProteinComplex','name'},[1,0]));
    % Enzyme number
    options.DnaA.number(iProtein) = iiProtein(1);
end

% Find segregation protein numbers in info structure
options.segregation.name = options.segregation.input(:,2);
options.segregation.type = options.segregation.input(:,1);
form = 'mature';
tProtein = length(options.segregation.name);
for iProtein = 1:tProtein
    % Find the Segregation enzyme numbers
    iiProtein = find(archive.infocmp(options.segregation.name{iProtein},{'info',options.segregation.type{iProtein},'name'},[1,0]));
    % Find the forms of the enzyme
    iiForm = {archive.info.(options.segregation.type{iProtein})(iiProtein).form};
    % Enzyme number of the correct form
    options.segregation.number(iProtein) = iiProtein(find(strcmpi(iiForm,form)));
end

%% ----------------- Processes state files ----------------- %

% ---- Get basic state files
state = archive.load_file(options.set,options.simulation,...
    'fileTypeTag',options.fileTypeTag);

% ---- Time
time = state.Time.values;
tTime = length(time);

% ---- DnaA
% Track the evolution of the all protein structures on chromosome
tmp = {state.Chromosome.complexBoundSites.vals};
clear data matrix
for iTime = 1:tTime;
    out = arrayfun(@(x)(sum(tmp{iTime}==x)),1:max(tmp{iTime}));
    matrix(1:length(out),iTime) = out;
end
% This matrix is a (pc by t) matrix where pc indicates the protein complex
% number and t the time point. the pc number corresponds with the pc number
% of the archive.info.ProteinComplex number.
state.processedChromosome.proteinComplexChromosome = matrix;
% Total number of proteins
tDnaA = length(options.DnaA.name);
% Which DnaA molecules should be included in the analysis
options.DnaA.useNumber = arrayfun(@(iDnaA)(max(strcmp(options.DnaA.name{iDnaA},{options.DnaA.useName{:}}))),1:tDnaA);
% Track the evolution of the DnaA complex from 2mers to 7mers
evolution = matrix(options.DnaA.number(options.DnaA.useNumber),:);
tUse = length(options.DnaA.useName);
% Check every time point the size of the largest DnaA complex
state.processedChromosome.evolveDnaAComplex = max((evolution>0) .* repmat([1:tUse],tTime,1)');

% ----- Segregation enzymes ----- %
clear matrix
tProtein = length(options.segregation.number);
for iProtein = 1:tProtein
    data = squeeze(state.(options.segregation.type{iProtein}).counts((options.segregation.number(iProtein)),1,:));
    matrix(iProtein,:) = data;
end
state.processedSegregation.matrix = matrix;
state.processedSegregation.evolution = sum(state.processedSegregation.matrix>1);

% ---- super coiling
% Get the number of coiled regions
state.processedChromosome.coiling = arrayfun(@(x)(size(state.Chromosome.superhelicalDensity(x).vals,1)),1:tTime);

%% ----------------- Determine regions ----------------- %

% Time and Volume
time = squeeze(state.Time.values);
volume = squeeze(state.Geometry.volume);

% ----- Pinching ----- %
% Color
analysis.pinching.color = 'r';
% y location
analysis.pinching.y = [0.4 1];
data = squeeze(state.Geometry.pinchedDiameter);
% Find the location of min an max value
maxValue=find(data==max(data));
minValue=find(data==0);
% If pinching is detected
if ~isempty(minValue)
    % Save time values
    analysis.pinching.time(1) = time(maxValue(end));
    analysis.pinching.time(2) = time(minValue(1));
    % Save min and max values
    analysis.pinching.value(1) = max(data);
    analysis.pinching.value(2) = min(data);
    % Set opaqueness
    analysis.pinching.opaqueness = options.alpha;
else
    analysis = rmfield(analysis,'pinching');
end

% ----- Gene duplication ----- %
% Color
analysis.gene_duplication.color = 'b';
% y location
analysis.gene_duplication.y = [0.4 1];
data = squeeze(sum(state.Chromosome.geneCopyNumbers,1));
% Find the location of min an max value
maxValue=find(data==max(data));
minValue=find(data==min(data));
% Save time values
if time(minValue(end)) < time(maxValue(1)); % If gene duplication acutally started
    analysis.gene_duplication.time(1) = time(minValue(end));
    analysis.gene_duplication.time(2) = time(maxValue(1));
    % Save min and max values
    analysis.gene_duplication.value(1) = min(data);
    analysis.gene_duplication.value(2) = max(data);
    % Set opaqueness
    analysis.gene_duplication.opaqueness = options.alpha;
else % If no gene duplication is detected
    % Remove field
    analysis = rmfield(analysis,'gene_duplication');
end


% ----- DnaA proteins ----- %
% Number of DnaA complex forms which are tracked
color = 'y';
opaquenessBase = 0.1;
y = [0,0.2];

analysis = get_regions(analysis,...
    state,'processedChromosome','evolveDnaAComplex',...
    'DnaA',color,opaquenessBase,y);

% ----- Segregation enzymes ----- %
color = 'r';
opaquenessBase = 0.1;
y = [0.2,0.4];

analysis = get_regions(analysis,...
    state,'processedSegregation','evolution',...
    'segregation_proteins',color,opaquenessBase,y);

% ----- Segregation completion ----- %
analysis.segregationCompletion.color = 'b';
analysis.segregationCompletion.opaqueness = options.alpha;
analysis.segregationCompletion.y = [0.9,1];
start = find(squeeze(state.Chromosome.segregated)==1,1,'first');
if ~isempty(start)
    analysis.segregationCompletion.time(1) = time(start);
    analysis.segregationCompletion.time(2) = time(find(squeeze(state.Chromosome.segregated)==1,1,'last'));
else
    analysis = rmfield(analysis,'segregationCompletion');
end


%% ----------------- Visualization processes ----------------- %

% New figure
if options.new
   figure; 
end

% Plot data
figObj = plot(time,volume);

% Get dimensions plot
options.xlim = xlim;
options.ylim = ylim;

% Plot regions
hold 'on'
field = fields(analysis);
tField = length(field);
for iField = 1:tField
    % Number of different sub regions
    tSubRegion = length(analysis.(field{iField}));
    for iSubRegion = 1:tSubRegion
    data = analysis.(field{iField})(iSubRegion).time;
    x = data([1,2,2,1]);
    height = analysis.(field{iField})(iSubRegion).y;
    y =  options.ylim([1,1,1,1]) + options.ylim(2)*[height(1),height(1),height(2),height(2)];
    color = analysis.(field{iField})(iSubRegion).color;
    opaqueness = analysis.(field{iField})(iSubRegion).opaqueness;
    p=patch(x,y,color,'EdgeColor','none');
    set(p,'FaceAlpha',opaqueness);
    end
end
hold 'off'

%% Set output arguments

varargout{1} = analysis;     % Contains the region information
varargout{2} = state;        % state file + additional information (processedX)

end

function analysis = get_regions(varargin)
% Add the different regions form a 1 by n array

analysis = varargin{1};
state = varargin{2};
field = varargin{3};
subfield = varargin{4};
fieldName = varargin{5};
color = varargin{6};
opaquenessBase = varargin{7};
y = varargin{8};

iSubRegion = 0;
time = squeeze(state.Time.values);

% Get max number of types of fields
tProtein = max(state.(field).(subfield));

for iProtein = 1:tProtein
    % See where the region is active
    tmp = find(state.(field).(subfield)==iProtein);
    if ~isempty(tmp)
        % Determine where the regions start
        tmpStart = ([1, (tmp(2:end) - tmp(1:end-1) ~= 1)]==1);
        tmpStart = tmp(tmpStart);
        % Determine where the regions stop
        tmpStop = [(tmp(1:end-1) - tmp(2:end) ~= -1),1]==1;
        tmpStop = tmp(tmpStop);
        % Number of regions found with this specific form of DnaA complex
        tRegion = length(tmpStart);
    else
        % If region does not exist
        tRegion = 0;
    end
    
    % Add the regions to analysis structure
    for iRegion = 1:tRegion
        iSubRegion = iSubRegion + 1;
        analysis.(fieldName)(iSubRegion).color = color;
        analysis.(fieldName)(iSubRegion).opaqueness = opaquenessBase*iProtein;
        analysis.(fieldName)(iSubRegion).y = y;
        try % Normal method
            timeStart = time(tmpStart(iRegion)-1);
            timeStop = time(tmpStop(iRegion));
        catch % Exception: This is required at t = 0
            timeStart = time(tmpStart(iRegion));
            timeStop = time(tmpStop(iRegion));
            if tmpStart(iRegion) == tmpStop(iRegion)
                timeStop = time(2)/2;
                analysis.(fieldName)(iSubRegion).opaqueness = 1;
            end
        end
        analysis.(fieldName)(iSubRegion).time(1) = timeStart;
        analysis.(fieldName)(iSubRegion).time(2) = timeStop;
    end
end

end