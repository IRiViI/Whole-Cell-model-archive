%Construct a detail stoichiometic structure from stoichiometic reaction
%formula in string form.
%
%   archive = extent_stoichiometry_archive(archive)
%
%   A stoichiometric formula is decomposed in its different components. The
%   compartments molecules, and reactions are saved as three different
%   fields. Where: The comparment field contains the id and name of the
%   different compartments. The molecule field contains the id, name,
%   its compartment, and the type of molecule (protein, metabolite etc).
%   The field reaction contains the number of the reaction according to the
%   "state" files of the simulations, id, name, reactionFormula (which is a
%   copy of the original stoichmetric reaction formula in string form), 
%   reversibility of the reaction and the molecules involved in the 
%   reaction.
%
%   NOTE: [1] Only molecules are included in the archive which are actually 
%   part of a reaction. [2] Molecules can be saved multiple times if they 
%   are active in different compartments.
%   
%   Generated fields: 
%   archive.stoichiometry
%   archive.info.MetabolicReaction.cytoscapID
%   
%   Required fields: 
%   archive.info.MetabolicReaction.Stoich  (reaction_stoichiometry_archive.m)
%
%   'progress'      - Show progress of function. 'On' or 'Off'. Default
%                     value is 'On'.

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 02/10/2016

function archive = extent_stoichiometry_archive(varargin)

% --------- Input Arguments --------- %

% Optional input arguments
if nargin > 1
    options = struct(varargin{2:end});
else
    options = struct();
end

% Process mandatory input arguments
archive = varargin{1};

% Add default settings if no values are given
if ~isfield(options,'process')
    options.process = 'On';
end

% --------- Settings --------- %

% Species types to add as nodes
speciesTypes = {'Metabolite','ProteinComplex','ProteinMonomer'};

% Compartments archive
compartmentList.id = {'c','d','e','m','tc','tm'};
compartmentList.name = {'Cytosol',...
    'DNA','Extracellular Space','Membrane','Terminal Organelle Cytosol'...
    'Terminal Organelle Membrane'};

% --------- Process --------- %

% Get stoichiometric matrix
archive = stoichiometric_matrix_maker(archive,compartmentList,speciesTypes,options);

end

function [archive] = stoichiometric_matrix_maker(varargin)
%Use the formula of the reaction to get an extensive structure about all
%   the reactions
%   archive.stoichiometry.molecule.(name,location,id,type)
%   archive.stoichiometry.reaction.molecule.(id,participation)

% Process mandatory input arguments
archive = varargin{1};
compartmentList = varargin{2};
speciesTypes = varargin{3};
options = varargin{4};

% Make structure for stoichiometric information
stoich = struct();
stoich.molecule = [];

% number of reactions
tReaction = length(archive.info.MetabolicReaction.Stoich);

% Progress spacing
if strcmp(options.process,'On')
    fprintf('\n     ');
end

% For every reaction
for iReaction = 1:tReaction
    
    % Show progress
    if strcmp(options.process,'On')
        display_progress(iReaction,tReaction);
    end
    
    % stoichiometic Formula
    stoichiometicFormula = archive.info.MetabolicReaction.Stoich{iReaction};
    
    % Get the individual elements of the reaction
    reactionElement = strsplit(stoichiometicFormula);
    
    % Check the reaction
    reactionStructure = check_elements(reactionElement,archive,speciesTypes);
    
    % ---------- Construct compartment list ----------  %
    
    % Add compartments
    for iCompartment = 1:length(compartmentList.id)
        stoich.compartment(iCompartment).id = compartmentList.id{iCompartment};
        stoich.compartment(iCompartment).name = compartmentList.name{iCompartment};
    end
    
    % ---------- Construct molecule list ----------  %
    
    % Two types of molecules
    moleculeField = {'Reactant','Product'};
    
    % For every molecule type
    for iField = 1:length(moleculeField)
        
        % For every molecule
        for iMolecule = 1:length(reactionStructure.(moleculeField{iField}).number);
            
            % Get the name, id and the compartment of the molecule
            nameMolecule = reactionStructure.(moleculeField{iField}).name{iMolecule};
            compartmentMolecule = reactionStructure.(moleculeField{iField}).location{iMolecule};
            idMolecule = reactionStructure.(moleculeField{iField}).id{iMolecule};
            type = reactionStructure.(moleculeField{iField}).type{iMolecule};
            tagMolecule = reactionStructure.(moleculeField{iField}).tag{iMolecule};
            
            % Check if the molecule is already used in a previously checked
            % reaction (checkMolecule = 0 or 1 for new or Pre-existing respectively)
            if iMolecule == 1 && iField == 1 && iReaction == 1
                checkMolecule = 0;      % Only for the first molecule
            else
                checkMoleculeID = strcmp({stoich.molecule(:).id},idMolecule);
                checkMoleculeComp = strcmp({stoich.molecule(:).compartment},compartmentMolecule);
                checkMolecule = sum(checkMoleculeID.*checkMoleculeComp);
            end
            
            if checkMolecule > 1
                error('Something is going wrong with molecule checking');
            end     
            
            % If molecule is new
            if checkMolecule == 0
                
                % Add to check list and stoich structure
                stoich.molecule(end+1).id = idMolecule;
                stoich.molecule(end).tag = tagMolecule;
                stoich.molecule(end).name = nameMolecule;
                stoich.molecule(end).compartment = compartmentMolecule;
                stoich.molecule(end).type = type;
            end
            
        end
        
    end
    
    % ---------- Construct reaction list ----------  %
    
    % Reaction number 
    stoich.reaction(iReaction).number = num2str(iReaction);
    
    % Reaction id
    reaction_id = archive.info.MetabolicReaction.ID{iReaction};
    stoich.reaction(iReaction).id = reaction_id;
    
    % Reaction name
    reaction_name = archive.info.MetabolicReaction.name{iReaction};
    % If there is no reaction name
    if isnan(reaction_name)
        reaction_name = '';
    end
    stoich.reaction(iReaction).name = reaction_name;
    
    % Reaction stoichiometic formula
    stoich.reaction(iReaction).reactionFormula = stoichiometicFormula;
    
    % Reaction reversibility
    stoich.reaction(iReaction).reversible = reactionStructure.reversible;
    
    % Reaction molecules (preallocation)
    stoich.reaction(iReaction).molecule = [];

    % Add reactants
    for iReactant = 1:length(reactionStructure.Reactant.number)
        stoich.reaction(iReaction).molecule(end+1).tag = reactionStructure.Reactant.tag{iReactant};
        stoich.reaction(iReaction).molecule(end).id = reactionStructure.Reactant.id{iReactant};
        stoich.reaction(iReaction).molecule(end).participation = 'reactant';
    end
    
    % Add products
    for iProduct = 1:length(reactionStructure.Product.number)
        stoich.reaction(iReaction).molecule(end+1).tag = reactionStructure.Product.tag{iProduct};
        stoich.reaction(iReaction).molecule(end).id = reactionStructure.Product.id{iProduct};
        stoich.reaction(iReaction).molecule(end).participation = 'product';
    end
    
    % Save reaction id list to the info list as well
    archive.info.MetabolicReaction.cytoscapeID{iReaction} = reaction_id;
    
end

% Save stoichiometic information
archive.stoichiometry = stoich;

end

function compartment = getLocation(element)
% Check if element indicates the location

switch element
    case '[c]:'
        compartment = 'c';
        case '[e]:'
        compartment = 'e';
    otherwise
        if strcmp(element(end),']') && strcmp(element(end-2),'[')
            compartment = 'different';
        else
            compartment = 0;
        end
end
end

function divider = getDivider(element)
% Check if element is a "+" sign or reaction vector "==>"

switch element
    case '+'
        divider = 'add';
    case '==>'
        divider = 'react';
    case '<==>'
        divider = 'revReact';    
    otherwise
        divider = 0;
end
end

function number = getNumber(element)
% Check if element is a number

if strcmp(element(1),'(') && strcmp(element(end),')')
    number = str2double(element(2:end-1));
    if number == 0
        error('One of the molecules in the reactions has a multiplier of zero in a reaction')
    end
else
    number = 0;
end
end

function  [molecule,type,name] = getMolecule(varargin)
% Check if element is a molecule.
%   The types of molecules to test are given in the speciesTypes list. The
%   lists them selves are given in the info structure. The element to
%   investigate is the element.

% Mandatory input arguments
element = varargin{1};
info = varargin{2};
speciesTypes = varargin{3};

% For every type of molecule
for iType = 1:length(speciesTypes)
    
    % Find for molecule
    moleculeNumber = find(strcmp(info.(speciesTypes{iType}).ID,element));
    
    % If molecule found
    if ~isempty(moleculeNumber)
        % Check if the molecule is found before
        if ~exist('molecule','var')
            % Save molecule number
            molecule = moleculeNumber;
            % Save the type of the molecule
            type = speciesTypes{iType};
            % Save the name of the molecule
            name = info.(type).name{moleculeNumber};
        else
            error('Molecule present in two different lists');
        end
    end
    
end

% If molecule does not exist
if ~exist('molecule','var')
    molecule = 0;
    type = 0;
    name = 0;
end

end

function reactionStructure = check_elements(reactionElement,archive,speciesTypes)
% Get a structure of the reactions that is easier to process

% Initiate structure
reactionStructure = struct();
reactionStructure.state = 'Reactant';

% counter for temperary number of product or reactant
PRcounter = 1;

% Get location
reactionStructure.compartment = getLocation(reactionElement{1});
if reactionStructure.compartment == 0
    reactionStructure.compartment = 'different';
%     errorMessage = sprintf('other location');
%     error(errorMessage)
end

for iElement = 1:length(reactionElement)
    
    % Check if element is location
    element = getLocation(reactionElement{iElement});

    if element(1) ~= 0 && ~strcmp(element,'different') % elements that contain both the location and the molecule are processed by the molecule indivication process
        element = 'location';
        if iElement ~= 1
            error('location other than first element')
        end
    else
        element = 0;
    end

    
    % Check if element is divider
    if element == 0
        element = getDivider(reactionElement{iElement});
        if element ~= 0
            switch element
                case 'add'
                    PRcounter = PRcounter + 1;
                case 'react'
                    reactionStructure.state = 'Product';
                    PRcounter = 1;
                    reactionStructure.reversible = 'false';
                  case 'revReact'
                    reactionStructure.state = 'Product';
                    PRcounter = 1;  
                    reactionStructure.reversible = 'true';
            end
        end
    end
    
    % Check if element is a number
    if element == 0
        element = getNumber(reactionElement{iElement});
        if element ~= 0
            reactionStructure.(reactionStructure.state).quantity(PRcounter) = element;
        end
    end
    
    % Check if element is a molecule
    if element == 0
        compartment = reactionStructure.compartment;
        if ~strcmp(compartment,'different')     % If the compartment is the same for all
            id = reactionElement{iElement};
            [element,type,name] = getMolecule(id,archive.info,speciesTypes);
        else                                    % If the compartment is different for different molecules
            id = reactionElement{iElement}(1:end-3); % Remove the location tag
            [element,type,name] = getMolecule(id,archive.info,speciesTypes);
            compartment = reactionElement{iElement}(end-1); % If the locations are not the same, there will be a tag at the end
            if ~(strcmp(compartment,'c') ||...
                    strcmp(compartment,'d') ||...
                    strcmp(compartment,'e') ||...
                    strcmp(compartment,'m') )
                % Check for some unexplected situation
               error('check compartment') 
            end
            
        end
        % If it's a molecule
        if element(1) ~= 0 
            
            % NOTE: Proteins have 5 different forms, only the first form is
            % used
            reactionStructure.(reactionStructure.state).number(PRcounter) = element(1);
            reactionStructure.(reactionStructure.state).name{PRcounter} = name;
            reactionStructure.(reactionStructure.state).id{PRcounter} = id;
            reactionStructure.(reactionStructure.state).tag{PRcounter} = [id '[' compartment ']'];
            reactionStructure.(reactionStructure.state).type{PRcounter} = type;
            reactionStructure.(reactionStructure.state).location{PRcounter} = compartment;
        end
    end
    
    % Exceptions (gamma_radiation is a stimulus)
    if strcmp(reactionElement{iElement}(1:end-3),'gamma_radiation')
        element = 'gamma_radiation';
        compartment = reactionElement{iElement}(end-1);
        name = 'gamma_radiation';
        reactionStructure.(reactionStructure.state).number(PRcounter) = 0;
        reactionStructure.(reactionStructure.state).name{PRcounter} = name;
        reactionStructure.(reactionStructure.state).id{PRcounter} = id;
        reactionStructure.(reactionStructure.state).tag{PRcounter} = [id '[' compartment ']'];
        reactionStructure.(reactionStructure.state).type{PRcounter} = 'radiation';
        reactionStructure.(reactionStructure.state).location{PRcounter} = compartment;
    end
    
    % If nothing is found
    if element == 0
        error('new type of element')
    end
    
end

end