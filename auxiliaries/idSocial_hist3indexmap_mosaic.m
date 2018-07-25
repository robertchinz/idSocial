function [nn,ctrs, edgesout,mapout,stdout,histarraytotal] = hist3indexmap_mosaic(varargin)
%[nn,ctrs,edgesout,mapout,stdout,histarray]=hist3indexmap_mosaic(ref_vector,'map',mval,'ctrs',ctrs,'histarrayedges',histarrayedges);
% hist3map is a modification of matlab's HIST3. In addition
% to the bivariate histogram, it is able to map the elements of input vector
% m which correspond to input vector x from which the histogram is produced
% into the corresponding bins.
% RCH 05/04/2010

% HIST3MAP(X)
%HIST3 Three-dimensional histogram of bivariate data.
%   HIST3(X) bins the elements of the M-by-2 matrix X into a 10-by-10 grid
%   of equally-spaced containers, and plots a histogram.  Each column of X
%   corresponds to one dimension in the bin grid.
%
%   HIST3(X,NBINS) plots a histogram using an NBINS(1)-by-NBINS(2) grid of
%   bins.  HIST3(X,'Nbins',NBINS) is equivalent to HIST3(X,NBINS).
%
%   HIST3(X,CTRS), where CTRS is a two-element cell array of numeric
%   vectors with monotonically non-decreasing values, uses a 2D grid of
%   bins centered on CTRS{1} in the first dimension and on CTRS{2} in the
%   second.  HIST3 assigns rows of X falling outside the range of that grid
%   to the bins along the outer edges of the grid, and ignores rows of X
%   containing NaNs.  HIST3(X,'Ctrs',CTRS) is equivalent to HIST3(X,CTRS).
%
%   HIST3(X,'Edges',EDGES), where EDGES is a two-element cell array
%   of numeric vectors with monotonically non-decreasing values, uses a 2D
%   grid of bins with edges at EDGES{1} in the first dimension and at
%   EDGES{2} in the second.  The (i,j)-th bin includes the value X(k,:) if
%
%      EDGES{1}(i) <= X(k,1) < EDGES{1}(i+1) and
%      EDGES{2}(j) <= X(k,2) < EDGES{2}(j+1).
%
%   Rows of X that that fall on the upper edges of the grid, EDGES{1}(end)
%   or EDGES{2}(end), are counted in the (I,j)-th or (i,J)-th bins, where
%   I and J are the lengths of EDGES{1} and EDGES{2}.  HIST3 does not count
%   rows of X falling outside the range of the grid.  Use -Inf and Inf in
%   EDGES to include all non-NaN values.
%
%   N = HIST3(X,...) returns a matrix containing the number of elements of
%   X that fall in each bin of the grid, and does not plot the histogram.
%
%   [N,C] = HIST3(X,...) returns the positions of the bin centers in a
%   1-by-2 cell array of numeric vectors, and does not plot the histogram.
%
%   HIST3(AX,X,...) plots into AX instead of GCA.
%
%   HIST3(..., 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   graphics parameter name/value pairs to fine-tune the plot.
%
%   Examples:
%
%      % Create the car data and make a histogram on a 7x7 grid of bins.
%      load carbig
%      X = [MPG,Weight];
%      hist3(X,[7 7]);
%      xlabel('MPG'); ylabel('Weight');
%
%      % Make a histogram with semi-transparent bars
%      hist3(X,[7 7],'FaceAlpha',.65);
%      xlabel('MPG'); ylabel('Weight');
%      set(gcf,'renderer','opengl');
%
%      % Make a histogram with bars colored according to height
%      hist3(X,[7 7]);
%      xlabel('MPG'); ylabel('Weight');
%      set(gcf,'renderer','opengl');
%      set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%
%      % Specify bin centers, different in each direction.  Get back
%      % counts, but don't make the plot.
%      cnt = hist3(X, {0:10:50 2000:500:5000});
%
%   See also ACCUMARRAY, BAR, BAR3, HIST, HISTC.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2009/05/07 18:31:14 $
tic
histarraytotal=[];
nn=[];
ctrs=[];
edgesout=[];
mapout=[];
stdout=[];
% spacing=7;

[cax,args,nargs] = axescheck(varargin{:});

if nargs < 1
    x=[1:10 10*ones(1,10) 1:10 ones(1,10) 5 6; 10*ones(1,10) 1:10 ones(1,10) 1:10 5 6]';
    %     x=[1:10 10*ones(1,9) 1:10 ones(1,9) 5 6; 9*ones(1,10) 1:9 ones(1,10) 1:9 5 6]';
    %     error('stats:hist3:TooFewInputs', 'Requires X.') % FOR TESTING;
    %     UNCOMMENT!!!
end
x = args{1};% FOR TESTING;
%     UNCOMMENT!!!

% try
[x sidx]=sortrows(x,1); % The plan is the following: First sort the Nx2 input
% catch 
%     keyboard
% end
orig_idx=sidx;
% vector by its first column, then produce maps for smaller ranges of
% values. Finally, put the maps together.
%
% if nargs > 1 && ~iscell(args{2})
%     m=args{2};
% else
%     m=[];
% end

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
pnames = {'nbins','ctrs','edges','map','histarrayedges','spacing'};
dflts =  { [],     [],       [],    [], [], 7};
% [errid,errmsg,nbins,ctrs,edges,map,histarrayedges, spacing,plotArgs] = internal.stats.getargs(pnames, dflts, args{:});
% if ~isempty(errmsg)
%     error(['stats:hist3:' errid], errmsg);
% end


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

%*disp(['[' mfilename '] ' ' Spacing of bins set to ' num2str(spacing)])
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
    map=map(sidx,:);
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


if ~isempty(nbins)
    % Use the specified number of bars in each direction, centers and edges
    % to be determined.
    histBehavior = true;
    if ~(isnumeric(nbins) && numel(nbins)==2)
        error('stats:hist3:BadNbins', ...
            'The number of bins must be specified with a 2-element numeric vector.');
    end
    autobins = true;
    
elseif ~isempty(ctrs)
    % Use the specified bin centers.
    histBehavior = true;
    if ~(iscell(ctrs) && numel(ctrs)==2 && isnumeric(ctrs{1}) && isnumeric(ctrs{2}))
        error('stats:hist3:BadCtrs', ...
            'Bin centers must be specified with a cell array containing two numeric vectors.');
    end
    ctrs = {ctrs{1}(:)' ctrs{2}(:)'};
    autobins = false;
    nbins = [length(ctrs{1}) length(ctrs{2})];
    
elseif ~isempty(edges)
    % Use the specified bin edges.
    histBehavior = false;
    if ~(iscell(edges) && numel(edges)==2 && isnumeric(edges{1}) && isnumeric(edges{2}))
        error('stats:hist3:BadEdges', ...
            'Bin edges must be specified with a cell array containing two numeric vectors.');
    end
    edges = {edges{1}(:)' edges{2}(:)'};
    autobins = false;
    % Just as with histc, there will be #edges bins
    nbins = [length(edges{1}) length(edges{2})];
    
else
    % Assume a 10x10 grid of bars, centers and edges to be determined.
    histBehavior = true;
    autobins = true;
    nbins = [20 20];
end

[nrows,ncols] = size(x);
if ncols ~= 2
    error('stats:hist3:WrongNumCols', 'X must be a matrix with two columns.');
end



% Special case for empty data (follows what HIST does).
if isempty(x)
    if autobins
        ctrs = {1:nbins(1) 1:nbins(2)};
    end
    n = zeros(nbins); % Nothing to count, return nbins(1) by nbins(2) zeros
    
else
    
    histcEdges=cell(2,1);
    binwidth=cell(2,1);
    xfilter1=true(1,size(x,1));
    xfilter2=true(1,size(x,1));
    for i = 1:2
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
            histcEdges{i} = [edges{i}(1)-binwidth{i}(1)*(1-1/spacing) edges{i}(1:end-1) edges{i}(end)+binwidth{i}(1)/spacing ]; % RCH: Changed edges{i}(2:end) to edges{i}(1:end); 2012/01/13
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
            histcEdges{i} = [edges{i}(1)-binwidth{i}(1)*(1-1/spacing) edges{i}(1:end-1) edges{i}(end)+binwidth{i}(1)/spacing ];
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
            %             histcEdges{i} = [edges{i}(1)-binwidth{i}(1)*(1-1/spacing) edges{i}(1:end-1) edges{i}(end)+binwidth{i}(1)*(1-1/spacing) ];
            histcEdges{i} = [edges{i}-binwidth{i}(1)*(1-1/spacing) edges{i}(end)+binwidth{i}(1)/spacing];% [edges{i}(1)-binwidth{i}(1)*(1-1/spacing) edges{i}(1:end-1) edges{i}(end)+binwidth{i}(1)*(1-1/spacing) ];
            
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
        
        x=x(xfilter1&xfilter2,:);
        orig_idx=orig_idx(xfilter1&xfilter2);
        if makemap; map=map(xfilter1&xfilter2,:); end;
        
    end
    
    %*disp(['[' mfilename '] ' ' Width of bins set to ' num2str([binwidth{1}(1) binwidth{2}(1)]) ' (x/y direction)'])
 try
    [av1, av2]=memory;
 catch
     av2.PhysicalMemory.Available = 8000000; %in bytes
 end
    maxmem=av2.PhysicalMemory.Available*.8;
    new_size=floor((maxmem-spacing*spacing*2*8)/(spacing*spacing*2*8)*(1/2));
 
    % Determine number of datapoints which fall in each bin
    size_cut=zeros(length(histcEdges{1})-1,1);
    for k=1:length(histcEdges{1})-1
        
        
        size_cut(k)=size(x(logical(x(:,1)>=histcEdges{1}(k)&x(:,1)<histcEdges{1}(k+1)),:),1);
        
    end
    if any((spacing*spacing*size_cut+spacing*spacing)*2*8>maxmem)
        keyboard;
        error('hist3map_mosaic:InsufficientMemory','Too much data in at least one bin.')
    end
    
    %     ntotal=zeros([ size(histcEdges{2},2)*spacing-spacing+1 size(histcEdges{1},2)*spacing-length(1:(cuts-1):length(histcEdges{1}))*spacing]);
    nsizetotal=[(spacing)*(size(histcEdges{2},2)-2)+spacing (spacing)*(size(histcEdges{1},2)-2)+spacing];
    ntotal=zeros(nsizetotal);
    
    for k=1:no_maps
        mapout{k}=NaN(size(ntotal));
        stdout{k}=NaN(size(ntotal));
        
    end
%     if makehists;
        histarraytotal=cell([size(ntotal,1) size(ntotal,2) 1]);
%     end;
    edgesout=cell(2,1);
    %     edgesout{1}=NaN(1,(size(histcEdges{1},2))*spacing-spacing);
    %     edgesout{2}=NaN(1,(size(histcEdges{2},2))*spacing-spacing);
    edgesout{1}=NaN(1,nsizetotal(2));
    edgesout{2}=NaN(1,nsizetotal(1));
    
    edgesouttemp=cell(2,1);
    mergeindex=1;
    %     for mosaic=1:(cuts-1):length(histcEdges{1})
    upbin=1;
    
    downbin=1;
    
    
    while upbin<length(size_cut)
        %     for mosaic=1:cuts+1
        %         sz = zeros(1,size(x,2));
        sz = zeros(1,size(x,2));
        
        histcEdgesMosaic=cell(2,1);
        %         histcEdgesMosaic{1}=histcEdges{1}(mosaic:1:min(mosaic+cuts,end));
        %         size(histcEdgesMosaic{2},2)*spacing size(histcEdgesMosaic{1},2)*spacing]
        %         if histcEdgesMosaic{1}(end)~=histcEdges{1}(end); histcEdgesMosaic{1}=[histcEdgesMosaic{1} histcEdges{1}(end)]; end;
        
        histcEdgesMosaic{2}=histcEdges{2};
        
        %         upbin=min(upbin+1,length(size_cut))
        
        while ~(sum(size_cut(downbin:min(upbin+1,length(size_cut))))>new_size) && ~(upbin+1>length(size_cut));
            upbin=min(upbin+1,length(size_cut));
            
        end
        
        %*disp(['[' mfilename '] ' 'Binrange ' num2str([downbin upbin]) ' of ' num2str(length(size_cut))] )
        
        histcEdgesMosaic{1}=histcEdges{1}(downbin:upbin+1);
        
        %         disp(['[' mfilename '] ' 'Limits ' num2str(histcEdgesMosaic{1}(1)) ':' num2str(histcEdgesMosaic{1}(end)) ' of ' num2str(histcEdges{1}(1)) ':' num2str(histcEdges{1}(end))] )
        
        
        xmosaic=x(logical(x(:,1)>=histcEdgesMosaic{1}(1)&x(:,1)<=histcEdgesMosaic{1}(end)),:);
        orig_idx_mosaic=orig_idx(logical(x(:,1)>=histcEdgesMosaic{1}(1)&x(:,1)<=histcEdgesMosaic{1}(end)));
        %         disp(['xm ' num2str(xmosaic(:,1)')])
        if makemap||makehists; mmosaic=map(logical(x(:,1)>=histcEdgesMosaic{1}(1)&x(:,1)<=histcEdgesMosaic{1}(end)),:); end;
        %         disp(['xm ' num2str(xmosaic(:,1)')])
        %         disp(['mm ' num2str(mmosaic(:,1)')])
        
        [mosaicnrows,mosaicncols]=size(xmosaic);
        %         bin = zeros(mosaicnrows,2);
        
        if ~isempty(xmosaic)
            bintotal=zeros(mosaicnrows*spacing,2);
            %     idbin=zeros(nrows*spacing,2);
            %             idx=perms(1:spacing);
            idx=1:spacing;
            %             permidx=[mod([0:size(idx,1)*size(idx,2)-1]',spacing)+1 reshape(idx',size(idx,1)*size(idx,2),1)];
            [p,q] = meshgrid(1:spacing, 1:spacing);
            permidx = [p(:) q(:)];
            %             psz=(spacing*factorial(spacing)*size(xmosaic,1)+spacing*factorial(spacing))*(2); %size(idx,1)*size(idx,2)*size(xmosaic,1);
            psz=(spacing*spacing*size(xmosaic,1)); %size(idx,1)*size(idx,2)*size(xmosaic,1);
            
            if (psz*2*8)>maxmem
                %*disp(['[' mfilename '] ' 'Size of index array: ' num2str(psz*2*8/(1024^2)) ' MB'])
                %*disp(['[' mfilename '] ' 'Maximum memory: ' num2str(maxmem/(1024^2)) ' MB'])
                keyboard;
                
                error('hist3map_mosaic:InsufficientMemory','Data split was not succesful, array would blow your computer.')
                
            else
                permidxtot=zeros(psz,2);
                permsize=whos('permidxtot');
                %*disp(['[' mfilename '] ' 'Size of index array: ' num2str(permsize.bytes/(1024^2)) ' MB'])
                %*disp(['[' mfilename '] ' 'Maximum memory: ' num2str(maxmem/(1024^2)) ' MB'])
            end
            spce=1;
            for k=1:size(xmosaic,1)
                %                 bintotal((spce-1)*size(permidx,1)+1:spce*size(permidx,1),:)=permidx(:,1:2)+spacing*(spce-1);
                permidxtot((spce-1)*size(permidx,1)+1:spce*size(permidx,1),:)=permidx(:,1:2)+spacing*(spce-1);
                spce=spce+1;
            end
            permidxtot=(permidxtot(permidxtot(:,1)>0,:));
            
            
            for i=1:2
                
                if i==1 ;edgesouttemp{i}=NaN(1,spacing*(upbin+1-downbin+1)); end;
                if i==2 ;edgesouttemp{i}=NaN(1,spacing*size(histcEdgesMosaic{i},2)); end;
                
                for spc=1:spacing
                    %                     disp([' upbin ' num2str(upbin) ' downbin ' num2str(downbin) ' i ' num2str(i) ' spc ' num2str(spc)]);
                    %                     edgesouttemp{i}(spc:spacing:end)
                    %                     histcEdgesMosaic{i}
                    histcEdgesTemp=histcEdgesMosaic{i}+(spc-1)*binwidth{i}(1)/spacing;
                    edgesouttemp{i}(spc:spacing:end)=histcEdgesTemp;
                    
                    
                    [~,bin] = histc(xmosaic(:,i),histcEdgesTemp,1);
                    %                     bin(bin<1,:)=2;
                    
                    sz(i)=length(histcEdgesTemp)-1;
                    
                    bin = (min(bin,length(histcEdges{i})-1));
                    
                    bintotal(spc:spacing:end,i)=(spacing)*(bin-1)+spc;% ((bin-1)*spacing+1+(spc-1))-1;
                    %                         if upbin==length(histcEdges{1})-1 && i==1; bintotal(spc:spacing:end,i)=((bin-1)*spacing+1+(spc-1)); keyboard;  end;
                    %                     if i==1
                    %                         disp(['Ed ' num2str(histcEdgesTemp)])
                    %                         disp(['bin ' num2str(bintotal(spc:spacing:end,i)')])
                    % %                         ((bin-1)*spacing+1+(spc-1))'
                    %
                    %
                    %                     end
                end % for spc=1:spacing
                
            end  %for i=1:2
            
            bintotal=([bintotal(permidxtot(:,2)',2) bintotal(permidxtot(:,1)',1)]);
            
            nsize=[(spacing)*(size(histcEdgesMosaic{2},2)-2)+spacing (spacing)*(size(histcEdgesMosaic{1},2)-1)+spacing];%[(size(histcEdgesMosaic{2},2))*spacing-spacing (size(histcEdgesMosaic{1},2))*spacing-spacing];
            bintotal=bintotal(all(permidxtot>0,2),:);
            %             bintotal=bintotal(all(bintotal'>0),:);
            
            n = accumarray(bintotal(all(bintotal'>0),:),1,nsize);
            
            %             clear bintotal;
            
            if makemap
                t=meshgrid(1:size(mmosaic,1),1:spacing*size(idx,2));
                t=reshape(t,1,size(t,1)*size(t,2))';
                orig_idx_mosaic=orig_idx_mosaic(t);
                for km=1:no_maps
                    
                    t2=mmosaic(t,km);
                    
                    
                    m{km} = accumarray(bintotal(all(bintotal'>0),:), t2(all(bintotal'>0),:),  nsize, @nanmean,NaN);
                    %                     m{km} = accumarray(bintotal, t2(all(bintotal>0, 2),:),  nsize,@min,NaN);
                    
                    
                    %                     st{km} = accumarray(bintotal, t2(all(bintotal>0, 2),:),  nsize, @nanstd,NaN);
                    st{km} = accumarray(bintotal(all(bintotal'>0),:), t2(all(bintotal'>0),:),  nsize, @nanstd,NaN);
                    % keyboard
                    
                end
            end
%             if makehists
%                 %                 histarray=NaN([size(n) histarraybins]);
%                 if ~makemap;
%                     t=meshgrid(1:length(mmosaic),1:spacing*size(idx,2));
%                     t=reshape(t,1,size(t,1)*size(t,2))';
%                     t2=mmosaic(t,km);
%                     
%                 end;
                bt=bintotal(all(bintotal'>0),:);
%                 t3=t2(all(bintotal'>0),:);
                orig_idx_mosaic=orig_idx_mosaic(all(bintotal'>0));
                
                %                 histarray = accumarray(bt,t3,[],@(x) {histc(x,histarrayedges)});
                histarray = accumarray(bt,orig_idx_mosaic,nsize,@(x) {x},{NaN}); % Found the idea for this at http://www.mathworks.com/matlabcentral/newsreader/view_thread/319998
%             end
            
            firstidx=1;
            ntotal(:,mergeindex:(mergeindex+nsize(2)-firstidx-spacing))=n(:,firstidx:end-spacing);
            
            
            %             edgesout{1}=0%(mergeindex:(mergeindex+nsize(2)-firstidx))=edgesouttemp{1}(firstidx:end-spacing);
            
            
            if makemap;
                for km=1:no_maps
                    mapout{km}(:,mergeindex:(mergeindex+nsize(2)-firstidx-spacing))=m{km}(:,firstidx:end-spacing);
                    stdout{km}(:,mergeindex:(mergeindex+nsize(2)-firstidx-spacing))=st{km}(:,firstidx:end-spacing);
                    
                end
                %                 figure
                %                 imagesc(mapout{1})
                %                 figure
                %                 imagesc(ntotal)
            end;
%             if makehists
                
                histarraytotal(:,mergeindex:(mergeindex+nsize(2)-firstidx-spacing))=histarray(:,firstidx:end-spacing);
%             end;
            ms=size(n(:,firstidx:end),2);
            %         mergeindex=mergeindex+size(histcEdgesMosaic{1},2)*spacing-spacing-1-spacing+1;
            mergeindex=mergeindex+ms-spacing;
            
            % Combine the two vectors of 1D bin counts into a grid of 2D bin
            % counts.
            %     n = accumarray(bin(all(bin>0,2),:),1,nbins);
        end
        if upbin==downbin
            upbin=upbin+1;
            downbin=downbin+1;
        else
            downbin=max(1,upbin);
            upbin=downbin+1;
        end
        
        
        
        
    end % while upbin~= etc..
    %     edgesout{2}=edgesouttemp{2}(1:end-spacing); % THE WHOLE EDGESOUT HISTORY IS BAAAAD ---WORK ON IT!!!
    %     edgesout{1}=edgesout{1}(floor(spacing/2)+1:end-floor(spacing/2));
    %     edgesout{2}=edgesout{2}(floor(spacing/2)+1:end-floor(spacing/2));
    edgesout=histcEdges;
    
    
    
    
    %*disp(['[' mfilename '] ' ' Needed ' num2str(toc) 's to calculate this.'])
    
    % keyboard
    if 0 < nargout
        %         nn = ntotal(spacing:end,spacing:end);
        nn = ntotal(spacing:end-(spacing),spacing:end-(spacing));
        %         nn=ntotal;
        nn=nn/(spacing*spacing); % Re-normalize,because now every data point can fall into spacing*spacing bins!
        if makemap
            for km=1:no_maps
                mapout{km,1}=mapout{km,1}(spacing:end-spacing,spacing:end-spacing);
                stdout{km,1}=stdout{km,1}(spacing:end-spacing,spacing:end-spacing);
                %                 mapout{km,1}=mapout{km,1}(1:size(nn,1),1:size(nn,2));
                %                 stdout{km,1}=stdout{km,1}(1:size(nn,1),1:size(nn,2));
            end
        end
%         if makehists
            
            histarraytotal=histarraytotal(spacing:end-spacing,spacing:end-spacing);
            
%         end
        return
    end
    
    if nargout == 0
        imagesc([edges{1}(1)-floor(spacing/2)*binwidth{1}(1)/spacing edges{1}(end)+floor(spacing/2)*binwidth{1}(1)/spacing],[edges{2}(1)-floor(spacing/2)*binwidth{2}(1)/spacing edges{2}(end)+floor(spacing/2)*binwidth{2}(1)/spacing],n);
        set(gca,'YDir','normal')
    end
end % if isempty(x)

