function [histog,bins]=histc_norm(datos,edges)

[histog,bins]=histc(datos,edges);
bw=diff(edges);
bw=bw(1);
bins=edges(1)+bw/2:bw:edges(end)-bw/2;
histog=histog/sum(histog);
histog=histog(1:end-1);



histog(1)=histog(1)/(bins(2)-bins(1));
histog(end)=histog(end)/(bins(end)-bins(end-1));
for c_bins=2:length(histog)-1
    histog(c_bins)=histog(c_bins)/((bins(c_bins+1)+bins(c_bins))/2 - (bins(c_bins)+bins(c_bins-1))/2);
end % c_bins
