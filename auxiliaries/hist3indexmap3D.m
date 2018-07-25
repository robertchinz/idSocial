function [nn,ctrs, edgesout,mapout,stdout,histarray] = hist3indexmap3D(varargin)
%[nn,ctrs,edgesout,mapout,stdout,histarray]=hist3indexmap_mosaic(ref_vector,'map',mval,'ctrs',ctrs,'histarrayedges',histarrayedges);
% hist3map is a modification of matlab's HIST3. In addition
% to the bivariate histogram, it is able to map the elements of input vector
% m which correspond to input vector x from which the histogram is produced
% into the corresponding bins.
% RCH 05/04/2010

tic
histarraytotal=[];
nn=[];
ctrs=[];
edgesout=[];
map=[];


[~,args,nargs] = axescheck(varargin{:});


x = args{1};
no_dim=size(x,2);

if no_dim<3
 
    edges{3}=[-inf inf];
end




% See if nbins/ctrs was given as the second argument, or only name/value
% pairs.
if nargs > 1 && ~ischar(args{2})
    binSpec = args{2};
    args = args(3:end); % strip off x and nbins/ctrs
else
    binSpec = [];
    args = args(2:end); % strip off x
end



% Process input parameter name/value pairs, assume unrecognized ones are
% graphics properties.
pnames = {'nbins','ctrs','edges','map','histarrayedges'};
dflts =  { [],     [],       [],    [], []};
ar=1;
while ar<=size(args,2)
    if ischar(args{ar}) && any(strcmp(args{ar},pnames))
        argtemp=args(ar+1);
        eval([args{ar} '=argtemp{1};'])
        pnames{strcmp(pnames,args{ar})}=[];
        ar=ar+2;
        
    else
        ar=ar+1;
    end
end

pnames_leftid=find(cellfun(@(x) ~isempty(x),pnames));
for pn=pnames_leftid
    eval([pnames{pn} '=[];'])
    pnames{pn}=dflts{pn};
end



% Make sure they haven't mixed 'nbins'/'ctrs'/'edges' name/value pairs with
% the CTRS or NBINS unnamed second arg syntax, or used more than one of
% those parameter name/value pairs.
if (isempty(nbins)+isempty(ctrs)+isempty(edges)+isempty(binSpec)) < 3
    error('stats:hist3:AmbiguousBinSpec', 'Ambiguous specification of bins.');
elseif ~isempty(binSpec)
    if iscell(binSpec)  % hist3(x,ctrs)
        ctrs = binSpec;
        %         keyboard
    else                % hist3(x,nbins)
        nbins = binSpec;
    end
end



if ~isempty(map)
    no_maps=size(map,2);
    makemap=true;
    if ~(isnumeric(map) || size(map,1)~=size(x,1))
        error('hist3map_mosaic:BadMapValues', ...
            'The mapping data has to have the same size as the primary data vector.');
    end
    m=cell(no_maps,1);
    st=cell(no_maps,1);
    mapout=cell(no_maps,1);
    stdout=cell(no_maps,1);
    
else
    makemap=false;
    mapout=[];
    stdout=[];
    no_maps=0;
end

if ~isempty(histarrayedges) && ~isempty(map) && nargout==6
    makehists=true;
    
elseif isempty(histarrayedges) && ~isempty(map) && nargout==6
    makehists=true;
    
    
else
    makehists=false;
    
    histarraytotal=[];
    histarray=[];
end


if ~isempty(edges)
    % Use the specified bin edges.
    histBehavior = false;
    if ~(iscell(edges) && numel(edges)>=2 && isnumeric(edges{1}) && isnumeric(edges{2}))
        error('stats:hist3:BadEdges', ...
            'Bin edges must be specified with a cell array containing two or three numeric vectors.');
    end
    edges = {edges{1}(:)' edges{2}(:)' edges{3}(:)'};
    autobins = false;
    % Just as with histc, there will be #edges bins
    nbins = [length(edges{1}) length(edges{2}) length(edges{3})];
    
else
    % Assume a 10x10 grid of bars, centers and edges to be determined.
    histBehavior = true;
    autobins = true;
    nbins = [20 20];
end

[nrows,ncols] = size(x);
if ncols < 2
    error('stats:hist3:WrongNumCols', 'X must be a matrix with two or three columns.');
end



% Special case for empty data (follows what HIST does).
if isempty(x)
    if autobins
        ctrs = {1:nbins(1) 1:nbins(2)};
    end
    n = zeros(nbins); % Nothing to count, return nbins(1) by nbins(2) zeros
    
else
    
    
    binwidth=cell(2,1);
    xfilter1=true(1,size(x,1));
    xfilter2=true(1,size(x,1));
    for i = 1:no_dim % if 3D, the z-coordinate will be added with a simple extra loop.
        % If only the number of bins was given, compute edges and centers
        % for equal-sized bins spanning the data.
        
        if autobins
            minx=min(x(:,i));
            maxx=max(x(:,i));
            if isinf(minx) || isinf(maxx)
                error('stats:hist3:InfData', ...
                    'Bin centers or edges must be specified when data contain infinite values.');
            elseif minx == maxx
                minx = minx - floor(nbins(i)/2) - 0.5;
                maxx = maxx + ceil(nbins(i)/2) - 0.5;
            end
            binwidth{i} = (maxx - minx) / (nbins(i)-1);
            edges{i} = minx + binwidth{i}*(0:nbins(i));
            %             edges{i}
            ctrs{i} = edges{i}(1:nbins(i)) + binwidth{i}/2;
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
        elseif histBehavior % if autobins
            c = ctrs{i};
            dc = diff(c);
            %             edges{i} = [c(1) c ] + [-dc(1) dc dc(end)]/2;
            edges{i} = [c(1) c c(end)+dc(end)] + [-dc(1) dc dc(end) dc(end)]/2;
            binwidth{i} = diff(edges{i});
            
            %             x(i,:)=x(x(:,i)>edges{i}(1),:);
            %             x(i,:)=x(x(:,i)<edges{i}(end),:);
            xfilter1=(x(:,i)>edges{i}(1));
            xfilter2=(x(:,i)<edges{i}(end));
            
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            % If the bin edges were given, compute their widths and centers (if
            % asked for).
        else % if ~histBehavior
            xfilter1=(x(:,i)>=edges{i}(1));
            xfilter2=(x(:,i)<edges{i}(end));
            e = edges{i};
            de = diff(e);
            binwidth{i} = diff(edges{i});
            % Make the display mimic bar's histc behavior: an implied bin
            % above edges(end), the same width as the last explicit one.
            % ctrs, edges, and binwidth need that explicitly, histcEdges
            % doesn't. RCH: SEEMS THAT IN MY CASE IT DOES?!!
            %             edges{i} = [e e(end)+de(end)];
            binwidth{i} = [de de(end)];
            
            if nargout > 1
                c = zeros(size(de));
                c(1) = e(1) + de(1)/2;
                for j = 2:length(c)
                    c(j) = 2*e(j) - c(j-1);
                end
                % When edges are specified, it may not be possible to return
                % centers for which the edges are midpoints.  Warn if that's
                % the case.
                if any(c <= e(1:end-1)) || ...
                        abs(c(end) - (e(end)-de(end)/2)) > 1000*eps(de(end));
                    warning('stats:hist3:InconsistentEdges', ...
                        'Cannot compute centers that are consistent with EDGES.');
                    c = e(1:end-1) + de/2;
                end
                ctrs{i} = [c e(end)+de(end)/2];
            end % if nargout > 1
        end % if autobins
        
        idx=1:nrows;
        x=x(xfilter1&xfilter2,:);
        orig_idx=idx(xfilter1&xfilter2);
        if makemap; map=map(xfilter1&xfilter2,:); end;
        
    end
    disp(['[' mfilename '] ' ' Width of bins set to ' num2str([binwidth{1}(1) binwidth{2}(1)]) ' (x/y direction)'])
    
    nrows_filtered=size(x,1);
    
    
    
    nsize=[(size(edges{2},2)-1) (size(edges{1},2)-1) (size(edges{3},2)-1)];
    ntotal=zeros(nsize);
    
    for k=1:no_maps
        mapout{k}=NaN(size(ntotal));
        stdout{k}=NaN(size(ntotal));
        
    end
    if makehists;
        histarraytotal=cell([size(ntotal,1) size(ntotal,2) 1]);
    end;
    
    
    bintotal=NaN(nrows_filtered,no_dim);
    for i=1:no_dim
        [~,bin] = histc(x(:,i),edges{i},1);
        bintotal(:,i) = (min(bin,length(edges{i})-1));
        
    end  %for i=1:dim
    
    
    
    n = accumarray(bintotal(all(bintotal'>0),:),1,nsize);
    
    %             clear bintotal;
    
    if makemap
        
        for km=1:no_maps
            
            t=map(:,km);
            m{km} = accumarray(bintotal(all(bintotal'>0),:), t(all(bintotal'>0),:),  nsize, @nanmean,NaN);
            st{km} = accumarray(bintotal(all(bintotal'>0),:), t(all(bintotal'>0),:),  nsize, @nanstd,NaN);
        end
    end
    if makehists
     
        
        histarray = accumarray(bintotal(all(bintotal'>0),:),orig_idx,nsize,@(x) {x},{NaN}); % Found the idea for this at http://www.mathworks.com/matlabcentral/newsreader/view_thread/319998
    end
    
    
    edgesout=edges;
    
    
    
    disp(['[' mfilename '] ' ' Needed ' num2str(toc) 's to calculate this.'])

    if 0 < nargout
        
        nn = n;
        if makemap
            for km=1:no_maps
                mapout{km,1}=m{km,1};
                stdout{km,1}=st{km,1};
                
            end
        end
        if makehists
            
            histarray=histarraytotal(1:size(nn,1),1:size(nn,2));
            
        end
        return
    end
    
    
end % if isempty(x)
