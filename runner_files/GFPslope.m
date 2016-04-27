classdef GFPslope < edu.stanford.covert.cell.sim.runners.SimulationRunner
    methods
        function this = GFPslope(varargin)
            this = this@edu.stanford.covert.cell.sim.runners.SimulationRunner(varargin{:});
        end
    end
    
    properties
        jobNumber = 0;              % Preallocate jobNumber
        addFeatureXmlFile = 'none';
        addFeatureXmlNumber = 'all';
    end
    
    methods (Access = protected)
        function modifyNetworkStructure(this, kb)
            
            % Import libraries
            import edu.stanford.covert.cell.sim.runners.functions.add_features;
            import edu.stanford.covert.cell.sim.runners.functions.get_info_kb;
            import edu.stanford.covert.cell.kb.Gene;
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.cell.kb.ProteinComplex;
            import edu.stanford.covert.cell.kb.Stimuli;
            import edu.stanford.covert.cell.kb.TranscriptionUnit;
            
            %             xmlDir = '/home/rick/Dropbox/Scripts/AdjustParamters/add_gene_GFP.xlsx';
            %             xmlDir = '/home/rickvink/WholeCell-master/src/+edu/+stanford/+covert/+cell/+sim/+runners/add_gene_GFP.xlsx';
            
            % Get mean half life
            meanHL = mean([kb.mRNAGenes.halfLife]);
            
            % Expression vector
            expression = zeros(1, 3);
            
            % Get Process_MacromolecularComplexation
            PMC = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');
            
            % Cytosol
            cytosol = kb.compartments(kb.cytosolCompartmentIndexs);
            
            % Extracellular space (atm not used)
            extracellularSpace = kb.compartments(kb.extracellularCompartmentIndexs);
            
            if ~strcmpi(this.addFeatureXmlFile,'none')
                % Get structure with new features
                [kb,~] = add_features(this.addFeatureXmlFile,...
                    'meanHL',meanHL,...
                    'expression',expression,...
                    'kb',kb,...
                    'PMC',PMC,...
                    'cytosol',cytosol,...
                    'setting_gene',this.addFeatureXmlNumber,...
                    'setting_proteinMonomer',this.addFeatureXmlNumber,...
                    'setting_proteinComplex',this.addFeatureXmlNumber,...
                    'setting_transcriptionUnit',this.addFeatureXmlNumber);
                
                % Make and save information structure
                info = get_info_kb(kb);
                [~,~,~] = mkdir(this.outDir);
                save([this.outDir '/' 'info'],'info');
                
            end
            
            % Set kb
            this.modifyNetworkStructure@edu.stanford.covert.cell.sim.runners.SimulationRunner(kb);
            
        end
        
        function modifyNetworkParameters(this, sim)
            
            % Change seed
            timeNow = clock;
            seed = this.jobNumber * 10000 + timeNow(2)*24*31 + timeNow(3)*24 + timeNow(4);  %Override seed
            this.seed = seed;
            
            % Get RNA information from xml file:
            
            % xml file
            %             afxmlf = '/home/rick/Dropbox/Scripts/AdjustParamters/add_gene_GFP.xlsx'; % Additional features xml file
            %             afxmlf = '/home/rickvink/WholeCell-master/src/+edu/+stanford/+covert/+cell/+sim/+runners/add_gene_GFP.xlsx';
            
            if ~strcmpi(this.addFeatureXmlFile,'none')
                
                % Get sheet names (fields)
                [~,fSheet] = xlsfinfo(this.addFeatureXmlFile);
                
                % Find RNA sheet
                iSheet = find(strcmpi(fSheet,'TranscriptionUnit'));
                
                % Get the raw data
                [~,~,raw] = xlsread(this.addFeatureXmlFile,iSheet);
                
                % Find the correct column
                cfelement = find(strcmpi(raw(1,:),'element'));
                iColumn = cfelement(3);
                
                % Number of rows
                tRow = size(raw,1);
                
                % Get list of rows that contains feature information
                tmp = {raw{:,1}}; % Get all the elements of the first row
                % Find al the elements that contain numerics
                premise1 = find(arrayfun(@(x)(isnumeric(tmp{x})),1:tRow));
                % Only select the selected transcription Units
                if ~isstr(this.addFeatureXmlNumber)
                    premise2 = premise1(find(arrayfun(@(x)(~isempty(find(tmp{x}==this.addFeatureXmlNumber))),premise1)));
                else
                    premise2 = premise1;
                end
                premise3 = premise2(premise2 > 0); % Empty/0 not allowed
                cfFeature = premise3;
                
                % select all
                select = 1:length(cfFeature);
                
                % List with all the selected RNA tags
                tagRNA = raw(cfFeature(select),iColumn);
                for i = 1:length(tagRNA)
                    tagRNA{i} = tagRNA{i}(2:end-1);
                end
                tagRNA = {tagRNA{:}};
                
                % edited original Add GFP script:
                g = sim.gene;
                time = sim.state('Time');
                rna = sim.state('Rna');
                trn = sim.process('Transcription');
                nascentMRNAIndexs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs, :), 1));
                [~, modTuIndexs] = ismember(tagRNA, rna.wholeCellModelIDs(rna.nascentIndexs));
                tuBindProb = trn.transcriptionUnitBindingProbabilities;
                meanBindProb = mean(tuBindProb(setdiff(nascentMRNAIndexs, modTuIndexs)));
                rnaDecayRates = rna.decayRates(rna.matureIndexs);
                modTuBindProb = tuBindProb;
                modTuBindProb(modTuIndexs) = meanBindProb;
                modTuBindProb(nascentMRNAIndexs) = modTuBindProb(nascentMRNAIndexs) * sum(tuBindProb(nascentMRNAIndexs)) / sum(modTuBindProb(nascentMRNAIndexs));
                modRnaExp = (rna.nascentRNAMatureRNAComposition * modTuBindProb) ./ (log(2) / time.cellCycleLength + rnaDecayRates);
                modRnaExp = modRnaExp / sum(modRnaExp);
                trn.transcriptionUnitBindingProbabilities = modTuBindProb;
                rna.expression(rna.matureIndexs) = modRnaExp;
            end
        end
    end
end