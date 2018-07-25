function idSocial_auxiliaries_setDynamicColormaps(fh,tickdensity_cm,colormap1,colormap2)
% http://www.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure

if nargin < 1 || isempty(fh)
    fh = gcf;
end
if nargin < 2 || isempty(tickdensity_cm)
    tickdensity_cm = 1;
end

if nargin < 3 || isempty(colormap1)
    load('colormap4InteractionMapsRedBlue')
    colormap1 = colormap4InteractionMapsRedBlue;
end
if nargin < 4 || isempty(colormap2)
    colormap2 = jet(64);
end

cm1_size = size(colormap1,1);

cm2_size = size(colormap2,1);

colormapComplete = [colormap2; colormap1];

img = findobj(fh,'Type','image','-not','Tag','TMW_COLORBAR');

if ~verLessThan('matlab','8.5.0')
    cbarhandle = findobj(fh,'Type','Colorbar');
else
    cbarhandle = findobj(fh,'Type','axes','Tag','Colorbar');
    
end


no_img = numel(img);
lims = NaN(no_img,2);
prts = NaN(1,no_img);
cdata = cell(no_img,1);
for k = 1:no_img
    cdata{k} = get(img(k),'CData');
    lims(k,1) = min(cdata{k}(:));
    lims(k,2) = max(cdata{k}(:));
    prts(k) = get(img(k),'Parent');
   
end
cm1 = find(lims(:,1)<0 & lims(:,2)>0);
cm2 = find(lims(:,1)>0 & lims(:,2)>0);



cm1lims = max([abs(min(lims(cm1,1))) max(lims(cm1,2))]);
cm1lims = [-cm1lims cm1lims];
cm2lims = [(min(lims(cm2,1))) max(lims(cm2,2))];
cm2lims = [0 max(lims(cm2,2))];

cdataScaled = cell(no_img,1);
for k = cm1'
    
    cdataScaled{k} = min(cm1_size,round((cm1_size-1)*(cdata{k}-cm1lims(1))/(cm1lims(2)-cm1lims(1)))+1);
    set(img(k),'CData',cdataScaled{k});
    
    caxis(prts(k),[1 size(colormapComplete,1)])
    
    set(cbarhandle(k),'Units','centimeter'); cbpos = get(cbarhandle(k),'Position');
    if strfind(lower(get(cbarhandle(k),'Location')),'north') || ...
            strfind(lower(get(cbarhandle(k),'Location')),'south')
        no_ticks = cbpos(3)*tickdensity_cm;
    else
        no_ticks = cbpos(4)*tickdensity_cm;
    end
    no_ticks = no_ticks/2+1;
    cbtickpos1 = 1:((cm1_size+1)/2-1)/(round(no_ticks-1)):(cm1_size+1)/2;
    cbtickpos2 = (cm1_size+1)/2:((cm1_size+1)/2-1)/(round(no_ticks-1)):(cm1_size);
    cbtickpos = [cbtickpos1 cbtickpos2(2:end)];
    
    
    % For labels:
    cbticklabel1 = cm1lims(1):(cm1lims(2))/(round(no_ticks-1)):0;
    cbticklabel2 = 0:(cm1lims(2))/(round(no_ticks-1)):cm1lims(2);
    cbticklabel = [cbticklabel1 cbticklabel2(2:end)];
    
    cblabelcell = cell(numel(cbtickpos),1);
    dtick = cbticklabel(2)-cbticklabel(1);
    
    dec =  find(num2str(dtick)=='.');
    prec1 = find(num2str(dtick)~='0' & num2str(dtick)~='.');
    prec = prec1(find(prec1>dec,1,'first'))-dec;
    if isempty(prec); prec=0; end;
    if any(prec1<dec) % Changes before decimal point
        cbtickpos = interp1(cbticklabel,cbtickpos,fix(cbticklabel));
        for lbc = 1:numel(cbtickpos)
            cblabelcell{lbc} = num2str(fix(cbticklabel(lbc)),'%.0f');
        end
    else
        for lbc = 1:numel(cbtickpos)
            cblabelcell{lbc} = num2str(cbticklabel(lbc),['%.' num2str(prec) 'f']);
        end
    end
    


    set(cbarhandle(k),'XTick',cbtickpos)
    if strfind(lower(get(cbarhandle(k),'Location')),'north') || ...
            strfind(lower(get(cbarhandle(k),'Location')),'south')
        set(cbarhandle(k),'XLim',[1 cm1_size])
        
        set(cbarhandle(k),'XTickLabel',cblabelcell)
    else
        set(cbarhandle(k),'YLim',[1 cm1_size])
        
    end
    
 end

for k = cm2'
    
    cdataScaled{k} = cm1_size + min(cm2_size,round((cm2_size-1)*(cdata{k}-cm2lims(1))/(cm2lims(2)-cm2lims(1)))+1);
    set(img(k),'CData',cdataScaled{k});
    caxis(prts(k),[1 size(colormapComplete,1)])
    
    set(cbarhandle(k),'Units','centimeter'); cbpos = get(cbarhandle(k),'Position');
    if strfind(lower(get(cbarhandle(k),'Location')),'north') || ...
            strfind(lower(get(cbarhandle(k),'Location')),'south')
        no_ticks = cbpos(3)*tickdensity_cm;
    else
        no_ticks = cbpos(4)*tickdensity_cm;
    end
    
    set(cbarhandle(k),'Units','centimeter'); cbpos = get(cbarhandle(k),'Position');
    cbtickpos = cm1_size+1:((cm2_size-1)/(round(no_ticks-1))):cm1_size+cm2_size;
    set(cbarhandle(k),'XTick',cbtickpos)
    
    % For labels:
    cbticklabel = cm2lims(1):(cm2lims(2)-cm2lims(1))/(round(no_ticks-1)):cm2lims(2);
    
    cblabelcell = cell(numel(cbtickpos),1);
    dtick = cbticklabel(3)-cbticklabel(2);
    dec =  find(num2str(dtick)=='.');
    prec1 = find(num2str(dtick)~='0' & num2str(dtick)~='.');
    prec = prec1(find(prec1>dec,1,'first'))-dec;
    if isempty(prec); prec=0; end;
    if any(prec1<dec) % Changes before decimal point
        cbtickpos = interp1(cbticklabel,cbtickpos,fix(cbticklabel));
        for lbc = 1:numel(cbtickpos)
            cblabelcell{lbc} = num2str(fix(cbticklabel(lbc)),'%.0f');
        end
    else
        for lbc = 1:numel(cbtickpos)
            cblabelcell{lbc} = num2str(cbticklabel(lbc),['%.' num2str(prec) 'f']);
        end
    end
    %     for lbc = 1:numel(cbtickpos)
%         dec =  find(num2str(dtick)=='.');
%         prec = find(num2str(dtick)~='0' & num2str(dtick)~='.',1,'first')-dec;
%         if isempty(prec); prec=0; end;
%         cblabelcell{lbc} = num2str(cbticklabel(lbc),['%.' num2str(prec) 'f']);
% 
%        
%     end
  
    if strfind(lower(get(cbarhandle(k),'Location')),'north') || ...
            strfind(lower(get(cbarhandle(k),'Location')),'south')
        set(cbarhandle(k),'XLim',[cm1_size+1 cm1_size+cm2_size])
        
        set(cbarhandle(k),'XTickLabel',cblabelcell)
    else
        set(cbarhandle(k),'YLim',[cm1_size+1 cm1_size+cm2_size])
        
    end
     

end


