function output_txt = myfunction(obj,event_obj,protein,PlottedProteinsIndex, Xorg, Yorg,PlottedRegionIndex,List0,scatterGroups,scatterGroupName)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
1;
% a=(scatterGroups==event_obj.Target); % protein types do not make different scatterGroups. Only Region(3) and phos (2). So maximum 6 groups. "a" tells to which scatterGroup does the current cursor belong to. 

s = get(event_obj,'Target');
scatterData = get(s, 'UserData');
PlottedProteinsIndex = scatterData.PlottedProteinsIndex;
PlottedRegionIndex = scatterData.PlottedRegionIndex;
Xorg = scatterData.Xorg;
Yorg = scatterData.Yorg;

dataIndex = get(event_obj,'DataIndex');
ProteinIndex=PlottedProteinsIndex(dataIndex);
SelectedRegion=List0{PlottedRegionIndex(dataIndex)};
% pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(Xorg(dataIndex))],...
    ['Y: ',num2str(Yorg(dataIndex))],protein(ProteinIndex).(SelectedRegion).LocusName,...
    [SelectedRegion ' length: ' num2str(protein(ProteinIndex).(SelectedRegion).RegionLength)]};


end
