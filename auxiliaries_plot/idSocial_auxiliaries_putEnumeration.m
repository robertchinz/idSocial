function th=idSocial_auxiliaries_putEnumeration(txt,enum_fontsize,ax,shift,alignment)

if nargin < 1 || isempty(txt)
    txt = 'EMPTY';
end

if nargin < 5 || isempty(alignment)
    alignment = {'bottom'; 'right'};
end

if nargin < 2 || isempty(enum_fontsize)
    enum_fontsize = 16;
end

if nargin < 3 || isempty(ax)
    ax = gca;
end

if nargin < 4 || isempty(shift)
    shift = [0 0];
end

orig_units = get(ax,'units');
set(ax,'Units', 'centimeter');
opos = get(ax,'Position'); tax = axes('units','centimeter',...
    'Position',[opos(1)+shift(1),opos(2)+opos(4)+shift(2), .5, .5]);
axis off
set(tax, 'XTick',[],'YTick',[], 'XTickLabel',[],'YTickLabel',[],'Visible','off')
th=text(0,0,.9,txt,'Fontsize',enum_fontsize,'VerticalAlignment',alignment{1},'HorizontalAlignment',alignment{2});
uistack(tax,'top')
set(ax,'Units', orig_units);
