% This function is mostly stolen from delay2grafo by Alfonso Pérez
% Escudero!
function [hf]=idSocial_displayGraph(delay,superanumbral,visible)

if nargin<2 || isempty(superanumbral)
    superanumbral=true(size(delay));
end
if nargin<3 || isempty(visible)
    visible='Off';
end

n_bichos=size(delay,1);

%%% Scaling
delay=delay/max(delay(:));

%%% Max. mean val makes top of hierarchy
delay_act=delay;
delay_act(~superanumbral)=NaN;
[~, pos]=sort(nanmean(delay_act,2),'ascend');
[~, pos]=sort(pos);

%%% Min. mean val makes top of hierarchy
[~, pos_min]=sort(nanmean(delay_act,2),'ascend');

x=0:1/(n_bichos-1):1;
% x=x.^.7;
hf=figure('Visible',visible);

colorines=[0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; .5 0 1; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 .5 .5; 0 1 1; 1 1 0;  .25 .5 .5; 1 .5 0; 0 .6 0; .5 .2 1; .5 1 1; 1 1 0];

radio=.08;
pos=pos-min(pos);
pos=pos/max(pos);
axis([-radio-.03 1+radio+.03 -radio-.03 1+radio+.03]);
hold on
for c1_bichos=1:n_bichos
    for c2_bichos=c1_bichos+1:n_bichos
        if superanumbral(c1_bichos,c2_bichos)
            if delay(c1_bichos,c2_bichos)>delay(c2_bichos,c1_bichos) % Arrow points to neighbour with negative lag with respect to focal. for < % ATTENTION: CHANGED DIRECTION! RCH
                inicio=[x(c1_bichos) pos(c1_bichos)];
                final=[x(c2_bichos) pos(c2_bichos)];
            else
                inicio=[x(c2_bichos) pos(c2_bichos)];
                final=[x(c1_bichos) pos(c1_bichos)];
            end
            vector=final-inicio;
            vector=vector/norm(vector);
            colorflechas=[.3 .3 .3];
            try
                %                 arrow_width=1.2-(abs(delay(c1_bichos,c2_bichos))-abs(posdelay_range(1)))/abs(diff(posdelay_range)); % Linewidth represents delay time: The longer the lag, the thinner the arrow
                %                 arrow_width=exp((posdelay_range(2)/delay(c1_bichos,c2_bichos))); % Linewidth represents delay time: The longer the lag, the thinner the arrow
                arrow_width=1;
                arrow(inicio+radio*vector,final-radio*vector,'LineWidth',3*arrow_width,'EdgeColor',colorflechas,'FaceColor',colorflechas,'TipAngle',10,'Length',20)
            catch
                keyboard
            end
        end % if supera umbral
    end % c2_bichos
end % c1_bichos

for c1_bichos=1:n_bichos
    try
        if all(isnan(pos)); pos=zeros(size(pos)); end;
        rectangle('Position',[x(c1_bichos)-radio pos(c1_bichos)-radio 2*radio 2*radio],'Curvature',[1 1],'LineWidth',3,'FaceColor',colorines(c1_bichos+1,:))
    catch
        keyboard
    end
    text(x(c1_bichos),pos(c1_bichos),num2str(c1_bichos),'FontSize',18,'FontWeight','bold','Color','k','HorizontalAlignment','center')
end % c1_bichos
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[],'XColor','w','YColor','w')
