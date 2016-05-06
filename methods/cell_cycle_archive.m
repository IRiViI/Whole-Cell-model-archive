% Example
%   cell_cycle_archive(archive,'set',1,'simulation',1)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 05/06/2016

function cell_cycle_archive(varargin)


% ----------------- Preparation processes ----------------- %

% Get archive
archive = varargin{1};

% Default settings
options.fileTypeTag = 'point_compressed_state.mat'; % State file name
options.alpha = 0.25;                               % Opaqueness of regions
options.new = true;                                 % New figure
options.info = 1;                                   % The information structure that should be used in this analysis
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

% Process inputs
% mandatory: 'set';'simulation'
tmpOptions = struct(varargin{2:end});
field = fields(tmpOptions);
for iField = 1:size(field,1)
    options.(field{iField}) = tmpOptions.(field{iField});
end

% Folder
folder = archive.set(options.set).simulation(options.simulation).folder;

% Check if the whole cell model tools are included
%--- Should still be made...

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

% Time
time = state.Time.values;

% Track the evolution of the all protein structures on chromosome
tmp = {state.Chromosome.complexBoundSites.vals};
tTime = length(time);
clear data
for iTime = 1:tTime;
    out = arrayfun(@(x)(sum(tmp{iTime}==x)),1:max(tmp{iTime}));
    matrix(1:length(out),iTime) = out;
end
% This matrix is a (pc by t) matrix where pc indicates the protein complex
% number and t the time point. the pc number corresponds with the pc number
% of the archive.info.ProteinComplex number.
state.processedChromosome.matrix = matrix;

% Total number of proteins
tDnaA = length(options.DnaA.name);
% Which DnaA molecules should be included in the analysis
options.DnaA.useNumber = arrayfun(@(iDnaA)(max(strcmp(options.DnaA.name{iDnaA},{options.DnaA.useName{:}}))),1:tDnaA);

% Track the evolution of the DnaA complex from 2mers to 7mers
evolution = matrix(options.DnaA.number(options.DnaA.useNumber),:);
tUse = length(options.DnaA.useName);
% Check every time point the size of the largest DnaA complex
state.processedChromosome.evolveDnaAComplex = max((evolution>0) .* repmat([1:tUse],tTime,1)');

% ----------------- Determine regions ----------------- %

% Time and Volume
time = squeeze(state.Time.values);
volume = squeeze(state.Geometry.volume);

% ----- Pinching ----- %
% Color
analysis.pinching.color = 'r';
% y location
analysis.pinching.y = [0.2 1];
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
analysis.gene_duplication.y = [0.2 1];
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


% % ----- DnaA proteins ----- %
% Number of DnaA complex forms which are tracked
tDnaA = max(state.processedChromosome.evolveDnaAComplex);
color = 'y';
opaquenessBase = 0.1;
iSubRegion = 0;
y = [0,0.2];
for iDnaA = 1:tDnaA
    % See where the region is active
    tmp = find(state.processedChromosome.evolveDnaAComplex==iDnaA);
    % Determine where the regions start
    tmpStart = ([1, (tmp(2:end) - tmp(1:end-1) ~= 1)]==1);
    tmpStart = tmp(tmpStart);
    % Determine where the regions stop
    tmpStop = [(tmp(1:end-1) - tmp(2:end) ~= -1),1]==1;
    tmpStop = tmp(tmpStop);
    % Number of regions found with this specific form of DnaA complex
    tRegion = length(tmpStart);
    % Add the regions to analysis structure
    for iRegion = 1:tRegion
        iSubRegion = iSubRegion + 1;
        analysis.dnaA(iSubRegion).color = color;
        analysis.dnaA(iSubRegion).opaqueness = opaquenessBase*iDnaA;
        analysis.dnaA(iSubRegion).y = y;
        try % Normal method
            timeStart = time(tmpStart(iRegion)-1);
            timeStop = time(tmpStop(iRegion));
        catch % Exception: This is required at t = 0
            timeStart = time(tmpStart(iRegion));
            timeStop = time(tmpStop(iRegion));
            if tmpStart(iRegion) == tmpStop(iRegion)
                timeStop = time(2)/2;
                analysis.dnaA(iSubRegion).opaqueness = 1;
            end
        end
        %         if tmpStart(iRegion) == tmpStop(iRegion)
        %             timeStart = timeStart - 1;
        %         end
        analysis.dnaA(iSubRegion).time(1) = timeStart;
        analysis.dnaA(iSubRegion).time(2) = timeStop;
    end
end
% % Find the location of min an max value
% maxValue=find(data==max(data));
% minValue=find(data==min(data));
% % % Find the location of min an max value
% maxValue=find(data==max(data,[],1));
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


end