%Add the stoichiometry formula to the archive.info.MetabolicReaction
%   structure. Additionally, information about enzymes and coenzymes are
%   added.
%
%   archive = stoichiometry_archive(archive,excelFile)
%
%   Stoichiometric reaction formula are obtained from the excel file. and
%   saved under "info" of the archive structure. The reactions are saved as a
%   single string.
%
%   Generated fields:
%   archive.info.MetabolicReaction.Stoich
%   archive.info.MetabolicReaction.Enzyme
%   archive.info.MetabolicReaction.Coenzyme
%
%   Required fields:
%   archive.info                   (labels_archive.m)
%
%   Required files:
%   Second input only (mmc4.xlsx)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 02/12/2016

% /home/rick/Documents/MIT_Project_Folder/Documents/Karr/mmc4.xlsx

function archive = stoichiometry_archive(varargin)

% ----- Parameters ----- %

% Process mandatory input argument
archive = varargin{1};
excelFile = varargin{2};

% Default settings
targets.StoichColumn = 9;               % Stoichiometry row of the sheet
targets.EnzymeIDColumn = 13;
targets.EnzymeCompartmentColumn = 14;
targets.Coenzyme = 15;
sheetTarget = 'S3O-Reactions';  % Name of sheet that contains the stoichimetry data

% ----- Process ----- %

% Get sheet names
[~,sheetName] = xlsfinfo(excelFile);

% Target sheet number
sheetNumber = find(strcmp(sheetName,sheetTarget));

% Get the raw data
[~,~,raw] = xlsread(excelFile,sheetName{sheetNumber});

% Get Stoichiometry data
archive = getStoichiometry(archive,raw,targets);

end

function archive = getStoichiometry(archive,raw,targets)
%Add the stoichiometry information out the raw data of excel sheet to archive.

% For every row of the sheet
for iRow = 1:size(raw,1)
    
    % Get the ID that is mentioned in the first column of that row
    ID = raw{iRow,1};
        
    % If the a reaction ID is found (or some other text)
    if ~isnan(ID)
        
        % Find reaction in metabolic reaction archive
        reaction = find(strcmp({archive.info.MetabolicReaction.ID},ID));
        
        % If reaction is found
        if ~isempty(reaction)
        
            % Add the requiested information.
            archive.info.MetabolicReaction(reaction).Stoich = raw{iRow,targets.StoichColumn};
            archive.info.MetabolicReaction(reaction).Enzyme.ID = raw{iRow,targets.EnzymeIDColumn};
            archive.info.MetabolicReaction(reaction).Enzyme.Compartment = raw{iRow,targets.EnzymeCompartmentColumn};
            archive.info.MetabolicReaction(reaction).Coenzyme = raw{iRow,targets.Coenzyme};
            
        end
    end    
        
end
end
