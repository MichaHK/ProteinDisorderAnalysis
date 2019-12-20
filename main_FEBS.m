function main_FEBS(handles)
% To Do: 
%The data cursor tip does not work if several regions are plotted together.
%Can be fixed by comparing the handle to "scatterGroup" with the propert
%event_obj.Target in the update function. I think it means that I need
%"Yorg" saved for each scatterGroup

persistent protein;

if (isempty(protein))
    %%
    %load('CustomTailData_20150425_223829.mat');
    load('CustomTailData_20150504_160000.mat'); % "Solvent accesssible" corrected
    
    hydrophobicityTable = FEBS_hydrophobicity(char([1:26] -1 + 'a'), [], 1);
    hydrophobicityTable = (hydrophobicityTable - 0.5); % Make -0.5 .. 0.5
    
    regions = {'head', 'rod', 'tail'};
    for i = 1:3
        for j = 1:numel(protein)
            protein(j).(regions{i}).CysNum = protein(j).(regions{i}).CysFraction * protein(j).(regions{i}).RegionLength;
            protein(j).(regions{i}).GlyNum = protein(j).(regions{i}).GlyFraction * protein(j).(regions{i}).RegionLength;
            protein(j).(regions{i}).CysPlusGlyFraction = protein(j).(regions{i}).CysFraction + protein(j).(regions{i}).GlyFraction;
            
            sequence = protein(j).(regions{i}).Sequence;
            sequenceNums = lower(sequence) - 'a' + 1;
            histogram = histc(sequenceNums, 1:26);
            protein(j).(regions{i}).HydrophobicityFractionPerAA = histogram .* hydrophobicityTable / protein(j).(regions{i}).RegionLength;
        end
    end
    
    1;
    if (0)
        %%
        which = [protein.IsKeratin] ~= 0;
        which = find(which);
        
        whichWhich = [cellfun(@(c)c(1)=='H', {protein(which).KeratinTissue})];
        which = which(whichWhich);
        
        hydrophobicityFractionPerAA = arrayfun(@(i)protein(i).tail.HydrophobicityFractionPerAA, which, 'UniformOutput', false);
        hydrophobicityFractionPerAA = vertcat(hydrophobicityFractionPerAA{:});
        
        figure(2); bar(hydrophobicityTable);
        set(gca, 'XTick', [1:26]);
        set(gca, 'XTickLabel', num2cell(char([1:26] - 1 + 'a')));
        title('hydrophobicity table');
        
        figure(3);
        bar(mean(hydrophobicityFractionPerAA, 1));
        set(gca, 'XTick', [1:26]);
        set(gca, 'XTickLabel', num2cell(char([1:26] - 1 + 'a')));
        title('hydrophobicity-fraction mean in HAIR keratin tails');
        
        figure(4);
        bar(std(hydrophobicityFractionPerAA, [], 1));
        set(gca, 'XTick', [1:26]);
        set(gca, 'XTickLabel', num2cell(char([1:26] - 1 + 'a')));
        title('std');
        
        %%
        figure(6);
        s = std(hydrophobicityFractionPerAA, [], 1);
        m = mean(hydrophobicityFractionPerAA, 1);
        candle(max(hydrophobicityFractionPerAA, [], 1)',...
            min(hydrophobicityFractionPerAA, [], 1)', (m + s)', (m - s)');
%         candle((m + 2*s)', (m - 2*s)', (m + s)', (m - s)');
        hold on;
        plot([0.5 26.5], [0 0], '-k');
        plot(m', '.k');
        hold off;
        set(gca, 'XTick', [1:26]);
        set(gca, 'XTickLabel', num2cell(char([1:26] - 1 + 'a')));
        1;
    end
end

try
    uiFieldNames = fieldnames(handles);
    whichGraphicHandles = cellfun(@(c)isscalar(handles.(c)), uiFieldNames);
    uiFieldNames = uiFieldNames(whichGraphicHandles);
    whichGraphicHandles = cellfun(@(c)ishghandle(handles.(c)), uiFieldNames);
    uiFieldNames = uiFieldNames(whichGraphicHandles);

    uiHandles = cellfun(@(c)double(handles.(c)), uiFieldNames);
    
    fieldsToSave = struct();
    fieldsToSave.edit = { 'String' };
    fieldsToSave.checkbox = { 'Value' };
    fieldsToSave.popup = { 'String', 'Value' };
    
%     data = struct();
%     fieldNamesToSave = fieldnames(fieldsToSave);
%     
%     GetFieldsValues = @(name, fields)cellfun(@(f)get(handles.(name), f), 'UniformOutput', false);
%     StoreFieldValues = @(name, fields, fieldValues)arrayfun(@(i)data.(name).(fields{i}) = fieldValues{i};
%     
%     StoreFieldValues(name, fieldsToSave.(uiType), GetFieldsValues(name, fieldsToSave.(uiType)))
    
catch
end





a=100;                                   % Surface of scatter points in figure.


List1={'Type', 'RegionLength', 'TailLength','HeadLength', 'FullLength', 'BodyLength', 'FractionPhosSites','PhosSitesPerPTMSites',...
        'FractionPTMSites','TailFracSize','HeadFracSize',...
        'JPredHelixPortion','JPredBetaPortion','JPredDisorderPortion', 'JPredBuried0Portion', 'JPredBuried5Portion','JPredBuried25Portion', ...
        'PsiHelixPortion','PsiBetaPortion','PsiCoilPortion', 'PsiDisorderScore', 'PsiDisorderPortion',...
        'ScratchSolAcc', 'ScratchSolAcc20', 'ScratchBuried', 'ScratchBuried20', 'SS3Helix', 'SS3Beta', 'SS3Other', 'SS8Helix', 'SS8Beta', 'SS8Other', ...
        'CysNum', 'CysFraction','GlyNum','GlyFraction', 'CysPlusGlyFraction', 'TailFraction','HeadFraction', 'BodyFraction'};  %Properties which do not depend on phosphorylation
List2={'TotalCharge','PartialCharge','PartialNegative','PartialPositive','RubinsteinChargeVar',...
    'ChargeSTD','AvgHydro','HydroSTD','HydroPhobicFrac','UverskyMeasure'}; %Properties that depend on phosphorylation
AllProperties=[List1 List2];

set(handles.PopUpX,'String',AllProperties);
set(handles.PopUpY,'String',AllProperties);
set(handles.filterPopup,'String',AllProperties);


typeToPlot=[];
for TypeInd=1:6             % 6 is the number of the known IF types
    checkBoxType=['CheckBoxType' '0'+TypeInd];
    typeToPlot=[typeToPlot get(handles.(checkBoxType),'Value')];
end


XSelectedIndFromAllProperties = get(handles.PopUpX, 'Value');
YSelectedIndFromAllProperties = get(handles.PopUpY, 'Value');

PropertyInList1Or2 = @(i)(1 + (i > numel(List1)));
XListNum = PropertyInList1Or2(XSelectedIndFromAllProperties);
YListNum = PropertyInList1Or2(YSelectedIndFromAllProperties);

XProperty = AllProperties{XSelectedIndFromAllProperties};
YProperty = AllProperties{YSelectedIndFromAllProperties};

%% Check Boxes: Regions to plot (head, rod, tail)
List0={'head' 'rod' 'tail'};
IsPlotHead=get(handles.HeadBox,'Value');
IsPlotRod=get(handles.RodBox,'Value');
IsPlotTail=get(handles.TailBox,'Value');

RegionsToPlot=[IsPlotHead IsPlotRod IsPlotTail];

%% Check Boxes: Type and Phosphorylation.
PhosStates={'Phos', 'NoPhos'};
IsPlotPhos=get(handles.PhosCheck,'Value');
IsPlotUnPhos=get(handles.UnPhosCheck,'Value');
if (XListNum==2 || YListNum==2)
    PhosStatesToPlot=[IsPlotPhos IsPlotUnPhos];
else
    PhosStatesToPlot=1;  %Irrelevant, just make that equals 1.
end

%%
if (get(handles.CheckPlotExternal,'Value'))
    figure(2)
    clf reset;
    hold on;
else
    axes(handles.axes1)
    cla reset
    hold on;
end

shouldDisplayLabels = get(handles.displayLabelsCheckbox, 'Value');

textLabels = [];
scatterGroups = [];
scatterGroupName=[];
x = [];
y = [];
c = [];
p = [];

shouldFilter = get(handles.filterCheckbox, 'Value');
filterPropertyIndex = get(handles.filterPopup, 'Value');
filteredProperties = get(handles.filterPopup, 'String');
filterProperty = filteredProperties{filterPropertyIndex};
filterRangeText = get(handles.filterRangeEditbox, 'String');
filter = @(x)1;

filterRange = eval(filterRangeText);

if (numel(filterRange) == 1)
    filter = @(x)x>filterRange;
elseif (numel(filterRange) == 2)
    filter = @(x)(x>=filterRange(1)) && x<=filterRange(2);
end



ShapeList={'s', '^', 'o'};
for RegionInd=1:numel(List0)
    SelectedRegion=List0{RegionInd};
    if RegionsToPlot(RegionInd)
        for PhosInd=1:numel(PhosStatesToPlot)
            if PhosStatesToPlot(PhosInd)==1
                PlotIndex=1;                    % we do not plot all proteins in every figure.
                for proteinIndex=1:numel(protein)
                    xCandidate=FixCandidate(protein,proteinIndex,XProperty,PhosStates,PhosInd,PhosStatesToPlot,XListNum,SelectedRegion);
                    xFullList(proteinIndex)=xCandidate;
                    yCandidate=FixCandidate(protein,proteinIndex,YProperty,PhosStates,PhosInd,PhosStatesToPlot,YListNum,SelectedRegion);
                    yFullList(proteinIndex)=yCandidate;
                    xOK=IsCandidateOK(xCandidate,XListNum);
                    yOK=IsCandidateOK(yCandidate,YListNum);
                    TypeOK=IsTypeOK(protein(proteinIndex).(SelectedRegion).Type,typeToPlot);

                    filterOK = 1;
                    if (shouldFilter)
                        filterCandidate = FixCandidate(protein,proteinIndex,filterProperty,PhosStates,PhosInd,PhosStatesToPlot,XListNum,SelectedRegion);
                        filterOK = filter(filterCandidate);
                    end
                    
                    if xOK*yOK*TypeOK*filterOK==1
                        x(PlotIndex) = xCandidate;
                        y(PlotIndex) = yCandidate;
                        c(PlotIndex)= protein(proteinIndex).(SelectedRegion).Type;
                        p(PlotIndex)= proteinIndex;
                        PlottedProteinsIndex(PlotIndex)=proteinIndex; % the index of the protein in the complete protein list.
                        PlottedRegionIndex(PlotIndex)=RegionInd;
                        PlotIndex=PlotIndex+1;
                    end
                end
                caxis([0,6])  % 6 is the number of the known IF types
                %         xAVG_all=mean(xFullList);
                %         xAVG_plotted=mean(x);
                %         xlabel([XProperty ' (#Proteins=' num2str(numel(x)) ')' ', AVG_{all}=' num2str(xAVG_all) ' ,AVG_{' num2str(numel(x)) '}=' num2str(xAVG_plotted) ])
                %         yAVG_all=mean(yFullList);
                xlabel(XProperty)

                %         yAVG_plotted=mean(y);
                %         ylabel([YProperty '  AVG_{all}=' num2str(yAVG_all) ' ,AVG_{' num2str(numel(y)) '}=' num2str(yAVG_plotted) ])
                ylabel(YProperty)
                CircleType={'filled', 'o'};
                %         DiamondType={'fill','o'};

                    Xorg=x;
                    Yorg=y;

                if get(handles.ABS_X_check,'Value')
                    x=abs(x);
                    X=x;
                    xlabel(['|' XProperty '|'])
                else
                    X=x;
                end
                if get(handles.ABS_Y_check,'Value')
                    y=abs(y);
                    Y=y;
                    ylabel(['|' YProperty '|'])
                else
                    Y=y;
                end

                %% Plot scatter
                s = scatter(X,Y,a,c,ShapeList{RegionInd},CircleType{PhosInd}, 'LineWidth',1.5); 
                scatterGroups(end+1) = s;
                
                scatterData = struct();
                scatterData.ProteinData = protein(p);
                
                scatterData.PlottedProteinsIndex = PlottedProteinsIndex;
                scatterData.PlottedRegionIndex = PlottedRegionIndex;
                scatterData.Xorg = Xorg;
                scatterData.Yorg = Yorg;
                
%                 scatterGroupName(end+1)={SelectedRegion PhosStates{PhosInd}};
                dcm_obj = datacursormode(handles.figure1);       
                set(dcm_obj,'UpdateFcn',{@myfunction, protein, PlottedProteinsIndex, Xorg, Yorg, PlottedRegionIndex, List0, scatterGroups, scatterGroupName})
            
                set(s, 'UserData', scatterData);
                
                %         yAid(1:numel(x))=yAVG_all;
                %         plot(x,yAid,'k')
                %         xAid(1:numel(y))=xAVG_all;
                %         plot(xAid,y,'k')

                %% Generate labels
                if (shouldDisplayLabels)
                    for pointIndex = 1:numel(X)
                        textLabels(end+1) = text(X(pointIndex), Y(pointIndex), ['  ' protein(p(pointIndex)).Gene], 'Rotation', 30, 'Interpreter', 'none', 'PickableParts', 'none', 'HitTest', 'off');
                    end
                end
            end

        end
    end
end
1;

%% Fix scales for graph
RoundTo = @(x, to)round(x ./ to) * to;
RoundConsistent = @(x)RoundTo(x, 10^(floor(log10(max(abs(x(:))))) - 2));

if get(handles.LogY_Check,'Value')
    y = [];
    for sg = 1:numel(scatterGroups)
        currentYValues = get(scatterGroups(sg), 'YData');
        y = [y; currentYValues(:)];
    end
    
    factorY = 1;
    if (min(abs(y(y~=0))) < 1)
        factorY = 1/min(abs(y(y~=0)))+1;
    end
    
    formula = @(x)sign(x).*log10(abs(x*factorY));
    
    for sg = 1:numel(scatterGroups)
        currentYValues = get(scatterGroups(sg), 'YData');
        set(scatterGroups(sg), 'YData', formula(currentYValues));
    end
    
    for i = 1:numel(textLabels)
        pos = get(textLabels(i), 'Position');
        pos(2) = formula(pos(2));
        set(textLabels(i), 'Position', pos);
    end
    
    YOrgTicks=get(gca,'ytick');
    YTicks=sign(YOrgTicks).*10.^(sign(YOrgTicks).*YOrgTicks)/factorY;
    YTicks=RoundConsistent(YTicks);
    set(gca,'YTickLabel',YTicks)
    
end

if get(handles.LogX_Check,'Value')
    x = [];
    for sg = 1:numel(scatterGroups)
        currentXValues = get(scatterGroups(sg), 'XData');
        x = [x; currentXValues(:)];
    end
    
    factorX = 1;
    if (min(abs(x(x~=0))) < 1)
        factorX = 1/min(abs(x(x~=0)))+1;
    end
    
    formula = @(x)sign(x).*log10(abs(x*factorX));
    
    for sg = 1:numel(scatterGroups)
        currentXValues = get(scatterGroups(sg), 'XData');
        set(scatterGroups(sg), 'XData', formula(currentXValues));
    end
    
    for i = 1:numel(textLabels)
        pos = get(textLabels(i), 'Position');
        pos(1) = formula(pos(1));
        set(textLabels(i), 'Position', pos);
    end
    
    XOrgTicks=get(gca,'xtick');
    XTicks=sign(XOrgTicks).*10.^(sign(XOrgTicks).*XOrgTicks)/factorX;
    XTicks=RoundConsistent(XTicks);
    set(gca,'XTickLabel',XTicks)
end

1;
if ~isempty(scatterGroups)
    LinePlotter(XProperty,YProperty,X,handles)
end

set(handles.TextPanel,'String', ['# plotted proteins: ' num2str(PlotIndex-1)])

%% Auxiliary functions
    function LinePlotter(XProperty,YProperty,x,handles)
        % axes(handles.axes1)
        % cla reset
        % hold on;
        if strcmp(XProperty,'AvgHydro')*strcmp(YProperty,'PartialCharge')==1
            if get(handles.LogY_Check,'Value')+get(handles.LogX_Check,'Value')==0
                if get(handles.ABS_Y_check,'Value')==1
                    q=2.785*x-1.151;
                    plot(x(q>0),q(q>0),'k')
                end
            end
        end
    end


    function TypeOK=IsTypeOK(Type,TypeToPlot)
        if TypeToPlot(Type)
            TypeOK=1;
        else
            TypeOK=0;
        end
    end

    function xCandidate=FixCandidate(protein,proteinIndex,XProperty,PhosStates,PhosInd,PhosStatesToPlot,XListNum,SelectedRegion)
        if numel(PhosStatesToPlot)==2
            if XListNum==2
                xCandidate=protein(proteinIndex).(SelectedRegion).(PhosStates{PhosInd}).(XProperty);
            else
                xCandidate=protein(proteinIndex).(SelectedRegion).(XProperty);
            end
        else
            xCandidate=protein(proteinIndex).(SelectedRegion).(XProperty);
        end
    end

    function xOK=IsCandidateOK(xCandidate,XListNum)
        if isempty(xCandidate)
            xOK=0; return
        end
        if XListNum==2
            xOK=1;
        end
        if XListNum==1
            if  xCandidate==0
                xOK=1;
            else xOK=1;
            end
        end
        
    end
end
