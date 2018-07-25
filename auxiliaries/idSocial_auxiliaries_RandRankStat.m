% function idSocial_auxiliaries_RandRankStat(n_delante,n_validos)
%% Permutations test to check the significance of the ranking for the three
%% days together
% I remove any fish that is disconnected in one of the days (it is only
% fish 4, the one that died shortly afterwards)
buenos=all(any(n_validos>100,2),3);
rangos=NaN(sum(buenos),3);
no_trials=size(n_delante,3);
for c_dias=1:no_trials
    medias=nanmean(n_delante(buenos,buenos,c_dias)./n_validos(buenos,buenos,c_dias),1);
    [s,orden]=sort(medias);
    rangos(orden,c_dias)=1:sum(buenos);
end
c=0;
difs=NaN(sum(buenos),3);
for c1_dias=1:3
    for c2_dias=c1_dias+1:3
        c=c+1;
        difs(:,c)=abs(rangos(:,c1_dias)-rangos(:,c2_dias));
    end % c2_dias
end % c1_dias
media_real=mean(difs(:));
% Ahora los casos random
medias_rand=NaN(1,100000);
for c_casos=1:100000
    rangos=[randperm(sum(buenos))' randperm(sum(buenos))' randperm(sum(buenos))'];
    c=0;
    difs=NaN(sum(buenos),3);
    for c1_dias=1:3
        for c2_dias=c1_dias+1:3
            c=c+1;
            difs(:,c)=abs(rangos(:,c1_dias)-rangos(:,c2_dias));
        end % c2_dias
    end % c1_dias
    medias_rand(c_casos)=mean(difs(:));
end % c_casos
figure
hist(medias_rand,20)
hold on
ejes=axis;
plot([media_real media_real],ejes(3:4),'k')
set(gca,'FontSize',15)
xlabel('Mean position difference')
ylabel('Number of cases')
p=sum(medias_rand<media_real)/length(medias_rand)
