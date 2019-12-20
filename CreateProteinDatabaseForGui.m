%% TO DO
% Deimination, also known as citrullination, is another
% form of post-translational modification of keratins and
% keratin filaments (Kubilus et al. 1980). During this process,
% the enzyme peptidylarginine deiminase converts arginine
% to citrulline. Because citrulline is neutral, the loss of a
% positive charge may lead to conformational changes of the
% keratin and keratin filaments

%% Data load
uniprotDataCacheFile = 'UniprotData.mat';
load('DataChart.mat');                                  % the 70 accession numbers are stored here.

if (~exist(uniprotDataCacheFile, 'file'))
    uniprotData = [];
    
    for i = 1:numel(dataChart.Entry)
        i
        a = getgenpept(dataChart.Entry{i});
        
        if (isempty(uniprotData))
            uniprotData = a;
        else
            uniprotData(i) = a;
        end
    end
    
    TypeList=importdata('TypeList.xlsx');       % Retrieve the IF type (1 to 6).
    for index=1:numel(uniprotData)
        name=uniprotData(index).LocusName;
        whichIs=find(cellfun(@(c) strcmp(name,c), TypeList.textdata(:,1)));
        uniprotData(index).Type=TypeList.data(whichIs);
    end
    
    save(uniprotDataCacheFile, 'uniprotData');
else
    load(uniprotDataCacheFile);
end

SelectedRegion=lower('tail.'); %always in lower case here.
% tailSequencesTable = {};
clear protein;
protein=struct;
%% Loop over Regions (Head, Rod, Tail)
RegionList={'head' 'rod' 'tail'}

for RegionInd=1:numel(RegionList)
    SelectedRegion=lower(RegionList{RegionInd});
    %% Loop over proteins.
    
    % for proteinIndex = 1:numel(uniprotData)
    for proteinIndex = 1:numel(uniprotData)
        %% Sequence for wanted region
        proteinIndex
        
        %a=getgenpept('Q9C075')
        a = uniprotData(proteinIndex);
        
        Regions = featuresparse(a,'feature','Region');
        notes = {Regions.note};                                     % list of regions (domains) in the whole protein
        whichIsRegionNote=find(cellfun(@(c)~isempty(strfind(lower(c), [SelectedRegion '.'])), notes));  % the "." prevents importing head1 and rod 2 etc'.
        
        gene = featuresparse(a,'feature','gene');
        geneName = gene.gene;
        protein(proteinIndex).Gene = gene.gene;
        
        if (numel(whichIsRegionNote) ~= 1)
            error('WRONG NUMBER OF TAILS!');
        end
        
        whichIsTailNote=find(cellfun(@(c)~isempty(strfind(lower(c), 'tail.')), notes));
        whichIsHeadNote=find(cellfun(@(c)~isempty(strfind(lower(c), 'head.')), notes));
        TailIndices=Regions(whichIsTailNote).Indices;
        HeadIndices=Regions(whichIsHeadNote).Indices;
        protein(proteinIndex).(SelectedRegion).TailFracSize=(TailIndices(2)-TailIndices(1)+1)/numel(a.Sequence);
        protein(proteinIndex).(SelectedRegion).HeadFracSize=(HeadIndices(2)-HeadIndices(1)+1)/numel(a.Sequence);
        
        
        
        
        RegionIndices=Regions(whichIsRegionNote).Indices;
        RegionSeq=a.Sequence(RegionIndices(1):RegionIndices(2));
        
        protein(proteinIndex).(SelectedRegion).LocusName=a.LocusName;
        protein(proteinIndex).(SelectedRegion).Region=SelectedRegion;
        protein(proteinIndex).(SelectedRegion).Accession=a.Accession;
        protein(proteinIndex).(SelectedRegion).Sequence=RegionSeq;
        protein(proteinIndex).(SelectedRegion).RegionLength=numel(RegionSeq);
        protein(proteinIndex).(SelectedRegion).TailLength=TailIndices(2)-TailIndices(1)+1;
        protein(proteinIndex).(SelectedRegion).HeadLength=HeadIndices(2)-HeadIndices(1)+1;
        protein(proteinIndex).(SelectedRegion).FullLength = numel(a.Sequence);
        protein(proteinIndex).(SelectedRegion).BodyLength = protein(proteinIndex).(SelectedRegion).FullLength ...
            - protein(proteinIndex).(SelectedRegion).HeadLength ...
            - protein(proteinIndex).(SelectedRegion).TailLength;
        protein(proteinIndex).(SelectedRegion).TailFraction = protein(proteinIndex).(SelectedRegion).TailLength / protein(proteinIndex).(SelectedRegion).FullLength;
        protein(proteinIndex).(SelectedRegion).HeadFraction = protein(proteinIndex).(SelectedRegion).HeadLength / protein(proteinIndex).(SelectedRegion).FullLength;
        protein(proteinIndex).(SelectedRegion).BodyFraction = protein(proteinIndex).(SelectedRegion).BodyLength / protein(proteinIndex).(SelectedRegion).FullLength;
        protein(proteinIndex).(SelectedRegion).RegionIndices=RegionIndices;
        protein(proteinIndex).(SelectedRegion).CysFraction=sum(lower(RegionSeq)=='c')/numel(RegionSeq);
        protein(proteinIndex).(SelectedRegion).GlyFraction=sum(lower(RegionSeq)=='g')/numel(RegionSeq);
        protein(proteinIndex).(SelectedRegion).CysPlusGlyFraction=protein(proteinIndex).(SelectedRegion).CysFraction + protein(proteinIndex).(SelectedRegion).GlyFraction;
        %% Parse Keratin data
        
        % To catch non- "KRT" prefix locus names that start with k
        if (a.LocusName(1) == 'K' && ~strcmp(geneName(1:3), 'KRT'))
            error('Bad "K" entry!');
        end
        
        protein(proteinIndex).IsKeratin = 0;
        
        % Is a Keratin?...
        if (numel(geneName) > 3 && strcmp(geneName(1:3), 'KRT'))
            % Schweizer, J., Bowden, P. E., Coulombe, P. a., Langbein, L., Lane, E. B., Magin, T. M., … Wright, M. W. (2006). New consensus nomenclature for mammalian keratins. Journal of Cell Biology, 174(2), 169–174. doi:10.1083/jcb.200603161
            
            protein(proteinIndex).IsKeratin = 1;
            
            if (geneName(end) >= 'A' && geneName(end) <= 'C')
                krtNum = str2double(geneName(4:end-1));
                krtNum = krtNum + 0.1 * (1 + (geneName(end) - 'A'));
            else
                krtNum = str2double(geneName(4:end));
            end
            
            
            1;
            protein(proteinIndex).KeratinNum = krtNum;
            if ((krtNum >= 1 && krtNum <= 28) || (krtNum >= 71 && krtNum <= 80))
                protein(proteinIndex).KeratinTissue = 'Epithelial';
                protein(proteinIndex).KeratinType = 1 + ~(krtNum >= 9 && krtNum <= 28);
            elseif ((krtNum >= 31 && krtNum <= 40) || (krtNum >= 81 && krtNum <= 86))
                protein(proteinIndex).KeratinTissue = 'Hair';
                protein(proteinIndex).KeratinType = 1 + (krtNum >= 81 && krtNum <= 86);
            else
                error('Bad Keratin entry!'); % Just in case
            end
        end
        
        
        %% Set type
        
        protein(proteinIndex).(SelectedRegion).Type=uniprotData(proteinIndex).Type;
        
        %% PTM sites
        Sites = featuresparse(a,'feature','Site');
        AllPhosList=[];
        TotalPTMinTail=0;
        SiteListInNotes={'Phosphotyrosine','Phosphoserine','Phosphothreonine','O-linked (GlcNAc)','N-acetylserine','N-acetylalanine'...
            'N6-acetyllysine','Cysteine methyl ester','Omega-N-methylarginine','Citrulline','N6-succinyllysine'};
        siteIdentifiers = cellfun(@(c)ToIdentifier(c), SiteListInNotes, 'UniformOutput', 0);
        
        if ~isempty(Sites)
            notes = {Sites.note};
            
            for i=1:numel(SiteListInNotes)
                % Check where they sit and whether they change the chrage....
                whichIs=find(cellfun(@(c)~isempty(strfind(c, SiteListInNotes{i})), notes));
                SiteAA_Indices=unique([Sites(whichIs).Indices]);
                SiteAA_Indices_in_Region=SiteAA_Indices(RegionIndices(1)<=SiteAA_Indices & SiteAA_Indices<=RegionIndices(2))-(RegionIndices(1)-1); % Indices match the region!!!!!!
                NumOfTailSites=length(SiteAA_Indices_in_Region);
                if ~isempty(strfind(lower(SiteListInNotes{i}),'phospho'))   % count phos sites for -2 charge.
                    AllPhosList=[AllPhosList SiteAA_Indices_in_Region];     % % Indices match the region, not the whole protein.
                end
                TotalPTMinTail=TotalPTMinTail+NumOfTailSites;
                
                protein(proteinIndex).(SelectedRegion).(['Sites__' siteIdentifiers{i}]) = SiteAA_Indices_in_Region;
                protein(proteinIndex).(SelectedRegion).(['Num' siteIdentifiers{i}]) = NumOfTailSites;
                protein(proteinIndex).(SelectedRegion).(['Fraction_' siteIdentifiers{i}]) = NumOfTailSites/protein(proteinIndex).(SelectedRegion).RegionLength;
                
            end
            protein(proteinIndex).(SelectedRegion).PhosSites = AllPhosList;
            protein(proteinIndex).(SelectedRegion).NumPhosSites = numel(AllPhosList);
            protein(proteinIndex).(SelectedRegion).FractionPhosSites = numel(AllPhosList)/protein(proteinIndex).(SelectedRegion).RegionLength;
            
            protein(proteinIndex).(SelectedRegion).NumPTMSites = TotalPTMinTail;
            protein(proteinIndex).(SelectedRegion).FractionPTMSites = TotalPTMinTail/protein(proteinIndex).(SelectedRegion).RegionLength;
            
            protein(proteinIndex).(SelectedRegion).PhosSitesPerPTMSites = numel(AllPhosList)/TotalPTMinTail;
        else
            protein(proteinIndex).(SelectedRegion).PhosSites = [];
            protein(proteinIndex).(SelectedRegion).NumPhosSites = 0; protein(proteinIndex).(SelectedRegion).FractionPhosSites=0;
            protein(proteinIndex).(SelectedRegion).NumPTMSites = 0;  protein(proteinIndex).(SelectedRegion).FractionPTMSites = 0;
            protein(proteinIndex).(SelectedRegion).PhosSitesPerPTMSites = 0;
            
            for i=1:numel(SiteListInNotes)
                protein(proteinIndex).(SelectedRegion).(['Num' siteIdentifiers{i}]) = 0;
                protein(proteinIndex).(SelectedRegion).(['Fraction_' siteIdentifiers{i}]) = 0;
            end
        end
        
        
        %% Charge
        PhosStates={'Phos','NoPhos'};
        for i=1:numel(PhosStates)
            if i==1
                PhosIndices=protein(proteinIndex).(SelectedRegion).PhosSites;
            else
                PhosIndices=[];
            end
            [protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AACharges, protein(proteinIndex).(SelectedRegion).(PhosStates{i}).TotalCharge, protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialCharge,...
                protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialNegative, protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialPositive, ...
                protein(proteinIndex).(SelectedRegion).(PhosStates{i}).ListNegSites, protein(proteinIndex).(SelectedRegion).(PhosStates{i}).ListPosSites]=FEBS_Charge(RegionSeq,PhosIndices);
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).ChargeSTD=std(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AACharges);
            %         DOI: 10.1051/jp2:1995157: (N+ - N-)^2 / (N+ + N-)  asymmetric
            %         when larger than 1, which promotes instability and expansion. See also fig 7 in
            %         10.1002/polb.20207. We changed to (TotalCharge)^2/(|NegCharge|+|PosCharge|)
            NumOfNegSites=sum(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).ListNegSites);
            NumOfPosSites=sum(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).ListPosSites);
            TotalCharge=protein(proteinIndex).(SelectedRegion).(PhosStates{i}).TotalCharge;
            TotalNeg=protein(proteinIndex).(SelectedRegion).RegionLength*protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialNegative;
            TotalPos=protein(proteinIndex).(SelectedRegion).RegionLength*protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialPositive;
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).RubinsteinChargeVar=sqrt((TotalCharge^2)/(abs(TotalNeg)+TotalPos));
            
            %         ((protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialCharge)^2)/...
            %             (protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialPositive+abs(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialNegative));
        end
        
        %% Hydrophobicity
        PhosStates={'Phos','NoPhos'};
        for i=1:numel(PhosStates)
            if i==1
                PhosIndices=AllPhosList;
            else
                PhosIndices=[];
            end
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AAHydros=FEBS_hydrophobicity(RegionSeq,PhosIndices);
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AvgHydro=sum(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AAHydros)/numel(RegionSeq);
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).HydroSTD=std(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AAHydros);
            NumOfHydrophobicAA=sum(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AAHydros>0.6);
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).HydroPhobicFrac=NumOfHydrophobicAA/protein(proteinIndex).(SelectedRegion).RegionLength;
            
            
        end
        
        %% Distance from Uversky 2D function (if >0 it is IDP).  PMID: 11025552
        % This is not loaded from uniprot data, but calculated based on previous KD and Charge
        for i=1:numel(PhosStates)
            protein(proteinIndex).(SelectedRegion).(PhosStates{i}).UverskyMeasure=abs(protein(proteinIndex).(SelectedRegion).(PhosStates{i}).PartialCharge)+1.151-2.785*protein(proteinIndex).(SelectedRegion).(PhosStates{i}).AvgHydro;
        end
        
        wmean = @(x, w)(w(:)'*x(:))/sum(w(:));
        
        %% Secondary structure + disorder prediction (from JPred)
        
        %load('JPredResults_20150406_101355.mat');
        load('JPredResults_20150414_141536.mat');
        proteinResult = ResultsByProtein.(protein(proteinIndex).(SelectedRegion).LocusName);
        whichTail = false(1, numel(proteinResult.Helix));
        whichTail(RegionIndices(1):RegionIndices(2)) = true;
        
        if (0)
            protein(proteinIndex).(SelectedRegion).JPredHelixPortion = mean(proteinResult.Helix(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBetaPortion = mean(proteinResult.Beta(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredDisorderPortion = mean(proteinResult.Disordered(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBuried0Portion = mean(proteinResult.Buried0(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBuried5Portion = mean(proteinResult.Buried5(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBuried25Portion = mean(proteinResult.Buried25(whichTail));
        else
            protein(proteinIndex).(SelectedRegion).JPredHelixPortion = wmean(proteinResult.Helix(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBetaPortion = wmean(proteinResult.Beta(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredDisorderPortion = wmean(proteinResult.Disordered(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBuried0Portion = wmean(proteinResult.Buried0(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBuried5Portion = wmean(proteinResult.Buried5(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).JPredBuried25Portion = wmean(proteinResult.Buried25(whichTail), proteinResult.Reliability(whichTail));
        end
        
        %% Secondary structure + disorder prediction (from PSIPRED)
        
        %load('PsiPredResults_20150424_225010.mat');
        load('PsiPredResults_20150429_101823.mat');
        proteinResult = ResultsByProtein.(protein(proteinIndex).(SelectedRegion).LocusName);
        whichTail = false(1, numel(proteinResult.Helix));
        whichTail(RegionIndices(1):RegionIndices(2)) = true;
        
        protein(proteinIndex).(SelectedRegion).PsiDisorderPortion = mean(proteinResult.Disordered(whichTail));
        protein(proteinIndex).(SelectedRegion).PsiDisorderPortionWeighed = wmean(proteinResult.Disordered(whichTail), proteinResult.Reliability(whichTail)); % Should not be used
        protein(proteinIndex).(SelectedRegion).PsiDisorderScore = mean(proteinResult.DisorderedScore(whichTail)); % Should not be used
        
        if (0)
            protein(proteinIndex).(SelectedRegion).PsiHelixPortion = mean(proteinResult.Helix(whichTail));
            protein(proteinIndex).(SelectedRegion).PsiBetaPortion = mean(proteinResult.Beta(whichTail));
            protein(proteinIndex).(SelectedRegion).PsiCoilPortion = mean(proteinResult.Coil(whichTail));
        else
            protein(proteinIndex).(SelectedRegion).PsiHelixPortion = wmean(proteinResult.Helix(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).PsiBetaPortion = wmean(proteinResult.Beta(whichTail), proteinResult.Reliability(whichTail));
            protein(proteinIndex).(SelectedRegion).PsiCoilPortion = wmean(proteinResult.Coil(whichTail), proteinResult.Reliability(whichTail));
        end
        
        %% Secondary structure + solvent accessibility (from SCRATCH)
        
        load('ScratchResults_20150425_222225.mat');
        proteinResult = ResultsByProtein.(protein(proteinIndex).(SelectedRegion).LocusName);
        whichTail = false(1, numel(proteinResult.SolvAcc));
        whichTail(RegionIndices(1):RegionIndices(2)) = true;
        
        protein(proteinIndex).(SelectedRegion).ScratchSolAcc = mean(proteinResult.SolvAcc20(whichTail) > 0.25);
        protein(proteinIndex).(SelectedRegion).ScratchSolAcc20 = mean(proteinResult.SolvAcc20(whichTail));
        
        protein(proteinIndex).(SelectedRegion).ScratchBuried = mean(proteinResult.SolvAcc20(whichTail) <= 0.25);
        protein(proteinIndex).(SelectedRegion).ScratchBuried20 = 1 - protein(proteinIndex).(SelectedRegion).ScratchSolAcc20;
        
        protein(proteinIndex).(SelectedRegion).SS3Helix = mean(proteinResult.SS3Helix(whichTail));
        protein(proteinIndex).(SelectedRegion).SS3Beta = mean(proteinResult.SS3Beta(whichTail));
        protein(proteinIndex).(SelectedRegion).SS3Other = mean(proteinResult.SS3Other(whichTail));
        
        protein(proteinIndex).(SelectedRegion).SS8Helix = mean(proteinResult.SS8Helix(whichTail));
        protein(proteinIndex).(SelectedRegion).SS8Beta = mean(proteinResult.SS8Beta(whichTail));
        protein(proteinIndex).(SelectedRegion).SS8Other = mean(proteinResult.SS8Other(whichTail));
        
    end
end

MatlabVersionInfo = ver('MATLAB');
if (str2double(MatlabVersionInfo.Version) >= 8.4)
    t = datetime('now');
    save(sprintf('CustomTailData_%i%02i%02i_%02i%02i%02i.mat', t.Year, t.Month, t.Day, t.Hour, t.Minute, round(t.Second)), 'protein');
else
    save('CustomTailData_20150504_160000.mat','protein')
end
