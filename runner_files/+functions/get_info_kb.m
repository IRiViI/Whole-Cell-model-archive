% Create an alternative info structure for the archive objects.
%
%   info = get_info_kb(kb)
%

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/24/2016

function out = get_info_kb(kb)

tMetabolites = length(kb.metabolites);
for iMetabolites = 1:tMetabolites
    out.Metabolite(iMetabolites).ID = kb.metabolites(iMetabolites).wholeCellModelID;
    out.Metabolite(iMetabolites).name = kb.metabolites(iMetabolites).name;
    out.Metabolite(iMetabolites).category = kb.metabolites(iMetabolites).category;
end

tMetabolicReaction = length(kb.reactions);
for iMetabolicReaction = 1:tMetabolicReaction
    out.MetabolicReaction(iMetabolicReaction).ID = kb.reactions(iMetabolicReaction).wholeCellModelID;
    out.MetabolicReaction(iMetabolicReaction).name = kb.reactions(iMetabolicReaction).name;
    out.MetabolicReaction(iMetabolicReaction).Stoich = nan;
    out.MetabolicReaction(iMetabolicReaction).Enzyme = nan;
    out.MetabolicReaction(iMetabolicReaction).Coenzyme = nan;
end

tProteinComplex = length(kb.proteinComplexs);
for iProteinComplex = 1:tProteinComplex
    out.ProteinComplex(iProteinComplex).ID = kb.proteinComplexs(iProteinComplex).wholeCellModelID;
    out.ProteinComplex(iProteinComplex).name = kb.proteinComplexs(iProteinComplex).name;
end

tProteinMonomer = length(kb.proteinMonomers);
for iProteinMonomer = 1:tProteinMonomer
    out.ProteinMonomer(iProteinMonomer).ID = kb.proteinMonomers(iProteinMonomer).wholeCellModelID;
    out.ProteinMonomer(iProteinMonomer).name = kb.proteinMonomers(iProteinMonomer).name;
end

tRna = length(kb.transcriptionUnits);
for iRna = 1:tRna
    out.Rna(iRna).ID = kb.transcriptionUnits(iRna).wholeCellModelID;
    out.Rna(iRna).name = kb.transcriptionUnits(iRna).name;
end

tStimulis = length(kb.stimulis);
for iStimulis = 1:tStimulis
    out.Stimulus(iStimulis).ID = kb.stimulis(iStimulis).wholeCellModelID;
    out.Stimulus(iStimulis).name = kb.stimulis(iStimulis).name;
end

end
