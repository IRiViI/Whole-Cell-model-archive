% Add new features to whole cell model simulations. It's a quite complex
%   proceedure. It's adviced to examine an example.
%
% [kb, features] = add_features('path/xml/file','kb',kb,...)
%
%   'setting_gene'      - Set which genes of the xml file should be added
%                         to genome and kb
%   'setting_transcriptionUnit'
%                       - Set which transcription units of the xml
%                         fileshould be added to genome and kb
%   'setting_proteinMonomer'
%                       - Similar as 'setting_gene' etc. for protein mono.
%   'setting_proteinComplex'
%                       - Similar as 'setting_gene' etc. for protein compl.
%
%   NOTE: The selection should be from 1 till a number x.
%   It's not possible to leave out a few features in between!
%   This because the area on the genome has to be specified on the genome.
%   It might actually be possible, but you should be careful!
%       Example: 'setting_gene',1:10 (possible) 'setting_gene',[3,5]
%       (probably results in problems)
%
%   Rules of the xml file:
%   - The sheet names should correspond with the name of the structures.
%   - The first column should contain the numbers.
%   - The first row should contain the labels of the row ('element')
%   - There may exist empty rows or rows with another purpose, as long as
%     these rows do not have a number in the first column.
%   - There may also be empty columns, as long as the the element of the
%     first row do not contain the label 'element'
%   - TRUE and FALSE must be presented as ISTRUE or ISFALSE in the xls file
%   - char string values must be given between ''. Example: {'protein'}
%     should be presented as 'protein' in the xml file
%   - Crosslinking within the xml file can be done using the *E:#* notation.
%     Where the E should refer to the structure and # to its number
%     Example: gene(1).transcriptionUnits = transcriptionUnit(5) should be
%     presented as: *transcriptionUnit:5* in the xml file.
%   - reference to variable without ''
%   - nan also without ''
%
%   Example (Simplified):
%    transcriptionUnit = TranscriptionUnit(...) % Left empty for this
%    example
%    gene = Gene(kb, NaN, {'Gene'}, 0, {''}, 580177, true, expression);
%    gene.transcriptionUnits = transcriptionUnit
%
%   XML file:
% field gene:
% number element element element element element element element element...
% 1      kb      nan    'Gene' 0       ''      580177 ISTRUE expression
%
% transcriptionUnits
% *transcriptionUnit:1*
%
% field transcriptionUnit (Left empty for the most part for this example!):
% number element ...
% 1

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/22/2016

function varargout = add_features(varargin)

% Process mandatory input argument
excelFile = varargin{1};

% Default settings
options.library =[];
options.setting_gene = 'all';
options.setting_proteinMonomer = 'all';
options.setting_proteinComplex = 'all';
options.setting_transcriptionUnit = 'all';

% Adjust options
inputOptions = struct(varargin{2:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

% Import libraries
import edu.stanford.covert.cell.kb.Gene;
import edu.stanford.covert.cell.kb.ProteinMonomer;
import edu.stanford.covert.cell.kb.ProteinComplex;
import edu.stanford.covert.cell.kb.Stimuli;
import edu.stanford.covert.cell.kb.TranscriptionUnit;

% Get sheet names (fields)
[~,fSheet] = xlsfinfo(excelFile);

% Number of sheets
tSheet = length(fSheet);

% Initiate xml file structure
parent = struct();

% Get genome
options.genome = options.kb.genome;

%--------------------------------------------
% Make element structures of the features and update genomic sequence
%--------------------------------------------

% Generate structure for every sheet For every sheet
for iSheet = 1:tSheet
    
    % name of current sheet
    nSheet = fSheet{iSheet};
    
    % Get the raw data
    [~,~,raw] = xlsread(excelFile,nSheet);
    
    % Create new structure
    parent.(nSheet) = struct();
    
    % Number of rows and columns
    tRow = size(raw,1);
    tColumn = size(raw,2);
    
    % Get labels
    fLabel = raw(1,1:tColumn);
    
    % Get list of column locations that contain the element label
    cfElement = find(strcmpi(fLabel,'element'));
    
    % Get list of rows that contains feature information
    tmp = {raw{:,1}}; % Get all the elements of the first row
    cfFeature = find(arrayfun(@(x)(isnumeric(tmp{x})),[1:tRow])); % Find al the elements that contain numerics
    cfFeature = cfFeature([tmp{cfFeature}]>0); % Empty/0 not allowed
    
    % Feature counter
    cFeature = 0;
    
    % Genome elements
    lGenome = find(arrayfun(@(x)(strcmp(fLabel{x}(1),'_')),[1:tColumn])); % column number
    cfGenome = fLabel(lGenome); % content
    
    % For every row that contains a feature
    for iRow = cfFeature
        
        % Add feature
        cFeature = cFeature + 1;
        
        % Element counter
        cElement = 0;
        
        % Clean Element
        clear('lElement');
        
        % For every column that is labeled with a gene label (ONLY check
        % for Gene)
        if length(cfGenome) > 0 && strcmpi(nSheet,'Gene')
            
            % Check if the gene should be included ( this is indeed quite
            % ugly.... but it works for now)
            % Look if the gene is in the given set.
            if ~isstr(options.setting_gene)
                tmp = find(options.setting_gene==raw{iRow,1});
                if isempty(tmp)
                    % If it's not found:
                    update_genome = 0;
                else
                    % If it's found:
                    update_genome = 1;
                end
            elseif strcmpi(options.setting_gene,'all')
                % Allways update if update_genome is all
                update_genome = 1;
            else
                error('Invalid input for setting_gene, only numerics and "''all''" is allowed')
            end
        else
            update_genome = 0;
        end
        
        if update_genome % If the genome should include this gene:
            try
                promoter = raw{iRow,lGenome(find(strcmp(cfGenome,'_genome:promoter')))}; % Get promoter
                promoter = promoter(2:end-1); % Remove the ''
            catch
                errorMessage = sprintf('No promoter found for gene %d in row %d of xml file',cFeature,iRow);
                error(errorMessage);
            end
            try
                gene = raw{iRow,lGenome(find(strcmp(cfGenome,'_genome:gene')))}; % Get promoter
                gene = gene(2:end-1); % Remove the ''
            catch
                errorMessage = sprintf('No gene found for gene %d in row %d of xml file',cFeature,iRow);
                error(errorMessage);
            end
            options.genome.sequence = [options.genome.sequence promoter gene];
            options.genome.sequenceLength = length(options.genome.sequence);
        end
        
        % For every column that contains an element
        for iColumn = cfElement
            
            % Count element
            cElement = cElement + 1;
            
            % Get element content
            element = raw{iRow,iColumn};
            
            % process element
            element = process_element(element,options);
            
            % Add element to element list
            lElement(cElement) = {element};
            
        end
        
        % Create feature
        switch nSheet
            case 'Gene'
                feature = Gene(lElement{:});
            case 'ProteinMonomer'
                feature = ProteinMonomer(lElement{:});
            case 'ProteinComplex'
                feature = ProteinComplex(lElement{:});
            case 'TranscriptionUnit'
                feature = TranscriptionUnit(lElement{:});
            case 'Stimuli'
                
        end
        
        % Add feature to parent
        try
            parent.(nSheet)(cFeature) = feature;
            parent.([nSheet '_info'])(cFeature).number = raw{iRow,1};
        catch
            parent.(nSheet) = feature;
            parent.([nSheet '_info']).number = raw{iRow,1};
        end
        
    end
    
end


%--------------------------------------------
% Process other
%--------------------------------------------

% Generate structure for every sheet For every sheet
for iSheet = 1:tSheet
    
    % name of current sheet
    nSheet = fSheet{iSheet};
    
    % Get the raw data
    [~,~,raw] = xlsread(excelFile,nSheet);
    
    % Number of rows and columns
    tRow = size(raw,1);
    tColumn = size(raw,2);
    
    % Get fields of labels
    fLabel = raw(1,1:tColumn);
    
    % Exception labels (labels starting with "_" are processed differently)
    lException = arrayfun(@(x)(strcmp(fLabel{x}(1),'_')),[1:tColumn]); % column number
    
    % Get the locations of the other labels
    cfLabel = find(~strcmpi(fLabel,'element') .*... % Exception: element
        ~strcmpi(fLabel,'label').*... % Exception: the name label
        ~lException); % Exception: a label starting with "_"
    
    % Feature counter
    cFeature = 0;
    
    % For every feature
    for iRow = cfFeature
        
        % Add feature
        cFeature = cFeature + 1;
        
        % Clean Element
        clear('lLabel');
        
        % For all the other labeled columns
        for iColumn = cfLabel
            
            % Get label name
            nLabel = fLabel{iColumn};
            
            % Get element content
            element = raw{iRow,iColumn};
            
            % Process element
            element = process_element(element,options,parent);
            
            % Add element
            parent.(nSheet)(cFeature).(nLabel) = element;
            
        end
    end
end


%--------------------------------------------
% Selection process
%--------------------------------------------

% Convert xml number to order number
lGene = get_selected_features(parent,options,'Gene','setting_gene');
% list that should be added to genome/kb structure
sGene = parent.Gene(lGene);
tGene = length(sGene);

% Convert xml number to order number
lTrans = get_selected_features(parent,options,'TranscriptionUnit','setting_transcriptionUnit');
% list that should be added to genome/kb structure
sTrans = parent.TranscriptionUnit(lTrans);
tTrans = length(sTrans);

% Convert xml number to order number
lProMo = get_selected_features(parent,options,'ProteinMonomer','setting_proteinMonomer');
% Add all desired protein monomers to genome
sProMo = parent.ProteinMonomer(lProMo);
tProMo = length(sProMo);

% Convert xml number to order number
lProCo = get_selected_features(parent,options,'ProteinComplex','setting_proteinComplex');
% Add all desired protein monomers to genome
sProCo = parent.ProteinComplex(lProCo);
tProCo = length(sProCo);

%--------------------------------------------
% Update genome and kb structure
%--------------------------------------------

genome = options.genome;

genome.genes(end+1:end+tGene) = sGene;
genome.transcriptionUnits(end+1:end+tTrans) = sTrans;

% Add genome to output structure (parent)
parent.genome = genome;

% Obtain kb
try
    kb = options.kb;
catch
    error(['kb structure is not given. please set kb.'...
        'Example: add_feature(''path/to/data'',...''kb'',kb,...'])
end

% Set features to kb structure
kb.genes = genome.genes;
kb.transcriptionUnits = genome.transcriptionUnits;
kb.proteinMonomers(end+1:end+tProMo) = sProMo;
kb.proteinComplexs(end+1:end+tProCo) = sProCo;


%--------------------------------------------
% Set output variables
%--------------------------------------------

varargout{1} = kb;
varargout{2} = parent;

end

function element = process_element(element,options,parent)
% Process Element

if ischar(element)
    if strcmpi(element,'nan')
        % Check if it is nan
        element = nan;
    elseif strcmpi(element,'istrue')
        % Check if it is true
        element = true;
    elseif strcmpi(element,'isfalse')
        % Check if it is false
        element = false;
    elseif strcmp(element,'''''')
        % If element is simply '' (order important)
        element = {''};
    elseif strcmp(element(1),'''') && strcmp(element(end),'''')
        % If the string is enclosed by '' (order important)
        element = {element(2:end-1)};
    elseif strcmp(element(1),'*') && strcmp(element(end),'*')
        % If the string is enclosed by **
        % Format: *field:number*
        % Example: *gene:4*
        %   element = parent.gene(4)
        tmp1 = element(2:end-1);
        tmp2 = strsplit(':',tmp1); % NOTE: This function is different than normal! (apparently...)
        % Find the correct field (regardless of captions and s at the end)
        fField = fields(parent);
        iField = find(strcmpi(fField,tmp2{1})); % Check if a match can be found regardless of caps
        if isempty(iField) % Check if a match can be found without the last letter (missing the "s")
            iField = find(strcmpi(fField,tmp2{1}(1:end-1)));
        end
        % Name of field
        nField = fField{iField};
        % Find correct number (look at the numbering of the xml sheet
        % (first column)
        iFeature = find([parent.([nField '_info']).number] == str2double(tmp2{2}));
        % Create new element
        element = parent.(nField)(iFeature);
    else
        % Else: It's a variable (order important)
        element = options.(element);
    end
end

end

function lGene = get_selected_features(varargin)
% Convert xml list number to order number

parent = varargin{1};
options = varargin{2};
type = varargin{3};
option = varargin{4};

% Total number of genes
xmlGene = [parent.([type '_info']).number];
tGene = length(parent.Gene);
if ~isstr(options.(option))
    % Add a selected few
    lGene = options.(option); % Get number set
    lGene = arrayfun(@(x)(sum(xmlGene(x)==lGene)),1:length(xmlGene)); % Check where the numbers of the given set are according to the xml file
    lGene = find(lGene); % Select the correct numbers
elseif strcmpi(options.(option),'all')
    lGene = 1:tGene;
else
    errorMessage = sprintf('Invalid input for %s, only numerics and "''all''" is allowed',option);
    error(errorMessage);
end

end