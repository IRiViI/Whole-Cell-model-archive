%Classify simulations in different groups. This version uses information
%   about the proteins and rnas to classefy the simulations.
%   archive = classification_archive(archive)
%
%   Generated fields:
%   archive.set.simulation.class
%   archive.set.class
%   archive.settings.classification
%
%   Required fields:
%   archive.set.simulation.mass
%   archive.set.simulation.lifetime
%   archive.set.simulation.trends
%
%   Required files:
%   None

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/05/2016

function archive = classification_archive(varargin)
% There are two parts to this function!
% First should the settings be specified.
% Then should this settings actually be implemented in the actual process


% Mandatory input argument
archive = varargin{1};

% Clear old settings
try
    archive.settings = rmfield(archive.settings,'classification');
catch
end

% -- Settings classification -- %

classification.Others = struct();

% First class: No growth. Mass at the end of the simulation is less than
% 120% as initial mass.
noGrowth.mass.max = 1.2;
classification.NoGrowth = noGrowth;

% Not used:

% % Second class: Stopped RNA production. Concentration of RNA goes down.
% stoppedRna.mass.min = 1.7;
% stoppedRna.rnaFactor.max = 1/50;
% classification.stoppedRna = stoppedRna;
%
% % Third class: Stopped protein production. Concentration of Proteins goes
% % down.
% stoppedProtein.mass.min = 1.7;
% stoppedProtein.rnaFactor.min = 1/50;
% stoppedProtein.proteinFactor.max = 1;
% classification.stoppedProtein = stoppedProtein;

% Class: Slow growth,
decreasingGrowth.lifetime.min = 50000;
decreasingGrowth.growth.max = 1.0;
classification.decreasingGrowth = decreasingGrowth;

% Class: No division while the mass is at least twice as high as at
% the beginning
% noDivision.mass.min = 1.7;
% noDivision.rnaFactor.min = 1/50;
% noDivision.proteinFactor.min = 1;
noDivision.lifetime.min = 50000;
noDivision.growth.min = 1.0;
classification.noDivision = noDivision;

% Class: Wild type like. It seems to behave quite reasonable
% vital.mass.min = 1.7;
% vital.rnaFactor.min = 1/50;
% vital.proteinFactor.min = 1;
vital.lifetime.max = 49999;
classification.vital = vital;

% Zeroth class: Others

% Add settings
archive.settings.classification = classification;


% -- Process classification -- %

% Total number of strains
tStrain = length(archive.set);

% Progress spacer
fprintf('Progress:\n      ')

% For every strain in the archive
for iStrain = 1:tStrain
    
    % Progres
    display_progress(iStrain,tStrain)
    
    % Total number of simulations
    tSimulation = length(archive.set(iStrain).simulation);
    
    % For every simulation of the strain
    for iSimulation = 1:tSimulation
        
        % Simulation parameters:
        mass1 = archive.set(iStrain).simulation(iSimulation).trends(1).mass;
        massEnd_3 = archive.set(iStrain).simulation(iSimulation).trends(end-3).mass;
        massEnd_2 = archive.set(iStrain).simulation(iSimulation).trends(end-2).mass;
        massEnd_1 = archive.set(iStrain).simulation(iSimulation).trends(end-1).mass;
        massEnd = archive.set(iStrain).simulation(iSimulation).trends(end).mass;
        %         rna1 = archive.set(iStrain).simulation(iSimulation).trends(1).rnas;
        %         rna3 = archive.set(iStrain).simulation(iSimulation).trends(3).rnas;
        %         protein1 = archive.set(iStrain).simulation(iSimulation).trends(1).proteins;
        %         proteinEnd = archive.set(iStrain).simulation(iSimulation).trends(end).proteins;
        timeEnd_3 = archive.set(iStrain).simulation(iSimulation).trends(end-3).time;
        timeEnd_2 = archive.set(iStrain).simulation(iSimulation).trends(end-2).time;
        timeEnd_1 = archive.set(iStrain).simulation(iSimulation).trends(end-1).time;
        timeEnd = archive.set(iStrain).simulation(iSimulation).trends(end).time;
        
        dmassEnd = (massEnd - massEnd_1)/(timeEnd - timeEnd_1);
        dmassEnd_1 = (massEnd_1 - massEnd_2)/(timeEnd_1 - timeEnd_2);
        dmassEnd_2 = (massEnd_2 - massEnd_3)/(timeEnd_2 - timeEnd_3);
        
        %         dmassEnd2 = (massEnd - massEnd_2)/(timeEnd - timeEnd_2);
        
        % Classification
        if  massEnd < noGrowth.mass.max * mass1
            
            % No growth
            class = 1;
            
            %         elseif massEnd >= stoppedRna.mass.min * mass1 &&...
            %                 rna3/rna1 < stoppedRna.rnaFactor.max
            %
            %             % Rna deficient
            %             class = 2;
            %
            %         elseif massEnd >= stoppedProtein.mass.min * mass1 &&...
            %                 rna3/rna1 > stoppedProtein.rnaFactor.min &&...
            %                 proteinEnd/protein1 < stoppedProtein.proteinFactor.max
            %
            %             % Protein deficient
            %             class = 3;
            %
        elseif ... %massEnd >= stoppedRna.mass.min * mass1 &&...
                (dmassEnd < decreasingGrowth.growth.max * dmassEnd_1 &&...
                dmassEnd_1 < decreasingGrowth.growth.max * dmassEnd_2) &&...
                timeEnd >= decreasingGrowth.lifetime.min
            
            % Slow growth
            class = 2;
            
        elseif ... %massEnd >= noDivision.mass.min * mass1 &&...
                ... %rna3/rna1 > noDivision.rnaFactor.min &&...
                ... %proteinEnd/protein1 > noDivision.proteinFactor.min &&...
                dmassEnd >= noDivision.growth.min * dmassEnd_1 &&...
                dmassEnd_1 >= noDivision.growth.min * dmassEnd_2 &&...
                timeEnd >= noDivision.lifetime.min
            
            % No division
            class = 3;
            
        elseif ... %massEnd >= vital.mass.min * mass1 &&...
                ... %rna3/rna1 > vital.rnaFactor.min &&...
                ... %proteinEnd/protein1 > vital.proteinFactor.min &&...
                timeEnd <= vital.lifetime.max
            
            % Vital
            class = 4;
            
        else
            
            % Other
            class = 0;
            
        end
        
        % Save class
        archive.set(iStrain).simulation(iSimulation).class = class;
        
    end
    
    % Get median class
    if ~isempty(archive.set(iStrain).simulation)
        subClassarchive = [archive.set(iStrain).simulation.class];
        classT = median(subClassarchive);
        if rem(classT,1) == 0
            archive.set(iStrain).class = classT;
        else
            classT = median(subClassarchive(subClassarchive ~= 0));
            if rem(classT,1) == 0
                archive.set(iStrain).class = classT;
            else
                warningMessage = sprintf('Set %d could not be well classified\n',iStrain);
                warning(warningMessage)
                archive.set(iStrain).class = 0;
            end
        end
    else
        warningMessage = sprintf('set %d does not contain any simulations',iStrain);
        warning(warningMessage)
    end
end

end
