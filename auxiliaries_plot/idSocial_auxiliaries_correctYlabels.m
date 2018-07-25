function idSocial_auxiliaries_correctYlabels(no_cols)
if nargin < 1 || isempty(no_cols)
    no_cols = 1;
end

haxes = findall(gcf,'type','axes');
h_axLabels = get(haxes,{'YLabel'});
stringTrue = cellfun(@(x) ~isempty(get(x,'String')),h_axLabels);
set([h_axLabels{:}],'Units','centimeter');
h_axLabels = h_axLabels(stringTrue,:);
h_parents = get([h_axLabels{:}] ,'Parent');
ylabelPos = cellfun(@(x) get(x,'Position'),h_axLabels,'UniformOutput',false );
ylabelPos = vertcat(ylabelPos{:});
ylabelPosX = ylabelPos(:,1);
ylabelPosY = ylabelPos(:,2);
ylabelPosZ = ylabelPos(:,3);

parentPos = cellfun(@(x) get(x,'Position'),h_parents,'UniformOutput',false );
parentPos = vertcat(parentPos{:});
parentPosX = parentPos(:,1);
parentPosY = parentPos(:,2);
parentPosZ = parentPos(:,3);



maxX = max(parentPosX);
minX = min(parentPosX);
% n = 2;
d = (maxX-minX)/no_cols;
if d==0 % All the same pos
    d=.1;
end
d2 = (maxX+d/2-(minX-d/2))/no_cols;
[ ~ , catYLabel] = histc(parentPosX,minX-d/2:d2:maxX+d/2);

for k = 1:no_cols
    idx = catYLabel==k;
    minLabelX = min(ylabelPosX(idx));
    for id = find(idx)'
        
%         disp([num2str(k) ',' num2str(id)])
%         get(h_axLabels{id},'String')
     set(h_axLabels{id},'Position',[minLabelX, ...
        ylabelPosY(id) ylabelPosZ(id)])
    end
end