function cell_cycle_archive(varargin)
% Example
%   cell_cycle_archive(archive,'set',1,'simulation',1)

% ----------------- Preparation processes ----------------- %

% Get archive
archive = varargin{1};

% Default settings
options.fileTypeTag = 'point_compressed_state.mat'; % State file name
options.alpha = 0.25;                               % Opaqueness of regions
options.new = true;                                 % New figure
options.info = 1;

options.DnaA.name = {'DnaA-ADP',... % Investigate these enzymes on the chromosome
    'DnaA-ATP',...
    'DnaA-ATP 2mer',...
    'DnaA-ATP 3mer',...
    'DnaA-ATP 4mer',...
    'DnaA-ATP 5mer',...
    'DnaA-ATP 6mer',...
    'DnaA-ATP 7mer',...
    };

% Process inputs
% mandatory: 'set';'simulation'
tmpOptions = struct(varargin{2:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

% Folder
folder = archive.set(options.set).simulation(options.simulation).folder;

% Find state files
lFile = find_files(folder,'state-');

% ----------------- Processes information ----------------- %

enzyme = options.DnaA.name;
tEnzyme = length(enzyme);
for iProtein = 1:tEnzyme
    % Find the DnaA enzymes
    iiProtein = find(archive.infocmp(enzyme{iProtein},{'info','ProteinComplex','name'},[1,0]));
    % Enzyme number
    options.DnaA.number(iProtein) = iiProtein(1);
end

% ----------------- Processes state files ----------------- %

% Get basic state files
state = archive.load_file(options.set,options.simulation,...
    'fileTypeTag',options.fileTypeTag);

% Add additional Chromosome information
tFile = length(lFile);

% For every state file
for iFile = 1:tFile
    
    % Load folder
    tmpState = load([folder '/' lFile(iFile).file],'Chromosome','Time');
    
    % ComplexBoundSites
    try
    [~,tmp] = find(tmpState.Chromosome.complexBoundSites(:,:,end));
    catch
    [~,tmp] = find(tmpState.Chromosome.complexBoundSites);
    end
    
    [a,b,~] = unique(tmp);
    
    % Find time point
    tmpTime = find(tmpState.Time.values(end)==state.Time.values);
        
    tEnzyme = length(options.DnaA.name);
    for iProtein = 1:tEnzyme
        iiProtein = options.DnaA.number(iProtein);
        tmp = find(a==iiProtein);
        if isempty(tmp)
            value = 0;
        else
            value = b(tmp);
        end
        state.extraChromosome.complexBoundSite(iProtein,1,tmpTime) = value;
    end
    
end

% ----------------- Main processes ----------------- %

% Time and Volume
time = squeeze(state.Time.values);
volume = squeeze(state.Geometry.volume);

% ----- Pinching ----- %
analysis.pinching.color = 'r';
data = squeeze(state.Geometry.pinchedDiameter);
% Find the location of min an max value
maxValue=find(data==max(data));
minValue=find(data==0);
% Warning
if isempty(minValue)
    warningMessage = sprintf('Pinching process did not stop according to data of set %d simulation %d',options.set,options.simulation);
    warning(warningMessage);
end
% Save time values
analysis.pinching.time(1) = time(maxValue(end));
analysis.pinching.time(2) = time(minValue(1));
% Save min and max values
analysis.pinching.value(1) = max(data);
analysis.pinching.value(2) = min(data);

% ----- Gene duplication ----- %
analysis.gene_duplication.color = 'b';
data = squeeze(sum(state.Chromosome.geneCopyNumbers,1));
% Find the location of min an max value
maxValue=find(data==max(data));
minValue=find(data==min(data));
% Save time values
analysis.gene_duplication.time(1) = time(minValue(end));
analysis.gene_duplication.time(2) = time(maxValue(1));
% Save min and max values
analysis.gene_duplication.value(1) = min(data);
analysis.gene_duplication.value(2) = max(data);

% % ----- Coiling ----- %
% analysis.gene_duplication.color = 'y';
% data = squeeze(sum(state.Chromosome.superhelicalDensity,1));
% % Find the location of min an max value
% maxValue=find(data==max(data));
% minValue=find(data==min(data));
% % Save time values
% analysis.gene_duplication.time(1) = time(minValue(end));
% analysis.gene_duplication.time(2) = time(maxValue(1));
% % Save min and max values
% analysis.gene_duplication.value(1) = min(data);
% analysis.gene_duplication.value(2) = max(data);

% ----------------- Visualization processes ----------------- %

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
    data = analysis.(field{iField}).time;
    x = data([1,2,2,1]);
    color =analysis.(field{iField}).color;
    p=patch(x,options.ylim([1,1,2,2]),color);
    set(p,'FaceAlpha',0.25);
end
hold 'off'


end