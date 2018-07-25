function idSocial_auxiliaries_correctXlabels(no_rows)

if nargin < 1 || isempty(no_rows)
    no_rows = 1;
end

haxes = findall(gcf,'type','axes');
h_axLabels = get(haxes,{'XLabel'});
stringTrue = cellfun(@(x) ~isempty(get(x,'String')),h_axLabels);
set([h_axLabels{:}],'Units','centimeter');
h_axLabels = h_axLabels(stringTrue,:);
h_parents = get([h_axLabels{:}] ,'Parent');
ylabelPos = cellfun(@(x) get(x,'Position'),h_axLabels,'UniformOutput',false );
ylabelPos = vertcat(ylabelPos{:});
ylabelPosX = ylabelPos(:,1);
ylabelPosY = ylabelPos(:,2);
ylabelPosZ = ylabelPos(:,3);

yExtent = get([h_axLabels{:}],'Extent');
yExtent = vertcat(yExtent{:});
yExtent = yExtent(:,2);

parentPos = cellfun(@(x) get(x,'Position'),h_parents,'UniformOutput',false );
parentPos = vertcat(parentPos{:});
parentPosX = parentPos(:,1);
parentPosY = parentPos(:,2);
parentPosZ = parentPos(:,3);



maxX = max(parentPosY);
minX = min(parentPosY);
d = (maxX-minX)/no_rows;
d2 = (maxX+d/2-(minX-d/2))/no_rows;
[ ~ , catYLabel] = histc(parentPosY,minX-d/2:d2:maxX+d/2);

for k = 1:no_rows
    idx = catYLabel==k;
    minLabelY = min(yExtent(idx));
    ExtentY = min(ylabelPosY(idx));
    for id = find(idx)'
        
%         disp([num2str(k) ',' num2str(id)])
     set(h_axLabels{id},'Position',[ylabelPosX(id) , ...
        minLabelY ylabelPosZ(id)],'VerticalAlignment','bottom')
    end
end