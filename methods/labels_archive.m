%Add the labels of metabolites, proteins, reactions etc. to the archive
%structure.
%
%   archive = labels_archive(archive,excelFile)
%
%   Alternative (if no archive file is yet created): 
%
%   archive = labels_archive(struct(),excelFile)
%
%   The labels for MetabolicReaction, Metabolite, ProteinComplex,
%   ProteinMonomer, Rna, Stimulus are extracted from the excel file
%   "excelFile". This is achieved by reading the different columns of the
%   sheets and saved under the titles "ID" and "name". Additonally are the
%   names of the different locations within the cell saved
%
%   Generated fields:
%   archive.info
%
%   Required fields:
%   None
%
%   Required files:
%   Second input only (StatePropertyRowColIDs.xlsx)

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 07/04/2016

% /home/rick/Documents/MIT_Project_Folder/Documents/Karr/StatePropertyRowColIDs.xlsx

function archive = labels_archive(varargin)


% Process mandatory input argument
archive = varargin{1};
excelFile = varargin{2};

% The locations of the statefiles
locations = {'Cytosol','Nucleoid','Extracellular Space','Membrane','Terminal Organelle cytosol','terminal Organelle membrane'};

% Get sheet names
[~,sheetName] = xlsfinfo(excelFile);

% For every sheet
for sheetNumber = 1: length(sheetName)
    
    % Get the raw data
    [~,~,raw] = xlsread(excelFile,sheetName{sheetNumber});
    
    % New clean structures
    ID = {};
    name = {};
    form = {};
    
    % Check if there is a type "form"
    if strcmpi(raw(1,4),'form')
        formState = true;
    else 
        formState = false;
    end
    
    % For every row of the table
    for rowNumber = 1:size(raw,1)
        
        % Get the number that is mentioned in that row
        number = raw{rowNumber,1};
        
        % Get the ID and name if the first element in the row is a number
        if isnumeric(number)
            id = raw(rowNumber,2);
            if ~isempty(id)
                if ~isnan(id{1})
                    ID(end+1) = id;
                    name(end+1) = raw(rowNumber,3);
                    % If there is also a form
                    if formState
                        form(end+1) = raw(rowNumber,4);
                    end
                end
            end
        end
        
    end
   
    % Get the first part of the sheet field name (required for saving)
    token = strtok(sheetName{sheetNumber},'.');
    
    % Save fields
    tEntry = length(ID);
    newField = struct();
    for iEntry = 1:tEntry
        newField(iEntry).ID = ID{iEntry};
        newField(iEntry).name = name{iEntry};
        if formState
            newField(iEntry).form = form{iEntry};
        end
    end
    archive.info(1).(token) = newField;
    
end

    % Add location info
    tLoc = length(locations);
    for iLoc = 1:tLoc
        archive.info.location(iLoc).name = locations{iLoc};
    end

end
