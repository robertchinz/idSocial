function trajectory_location = idSocial_auxiliaries_createTrajectoryLocation(directory,fname,parse_and_order)

if nargin<3 || isempty(parse_and_order)
    parse_and_order=[];
end

trajectory_location = idSocial_recursiveFolderSearch(directory,fname);
trajectory_location = idSocial_recursiveCountElemsInBranches(trajectory_location);
[elist_all, pathlist]= idSocial_GetCellDepth(trajectory_location);
% trajectory_location = idSocial_recursiveCutEmptyBranches(trajectory_location);

depthlist=cellfun(@(x) size(x,2),elist_all);
maxdepth=max(depthlist);
elist = NaN(size(elist_all,1),maxdepth);
for k=1:size(elist_all,1)
    elist(k,1:size(elist_all{k},2))=elist_all{k};
end

% all(repmat(elist(1,:),[size(elist,1),1])==elist,1)
superflous_subdivision=find(all(1==elist,1),1);


%% Get dimension names (path names)


shortpathlist=cellfun(@(x) strrep(x,directory,''),pathlist,'UniformOutput',false);
parsepath=cellfun(@(x) textscan(x,'%s','Delimiter',['\'  filesep]),shortpathlist,'UniformOutput',false);
parsepath=vertcat(parsepath{:});
parsepath= cellfun(@(x) x',parsepath,'UniformOutput',false);
parsepath=vertcat(parsepath{:});
parsepath=parsepath(:,1:end-1); % remove 'trajectories.mat'-string


parsepath(:,superflous_subdivision+1)=[];
dimnames=cell(1,size(parsepath,2));
for col=1:size(parsepath,2)
    dimnames{col}=unique(parsepath(:,col));
end

dimspecs = cell(size(dimnames,2),2);
numerated_charspec= NaN(size(dimnames,2),1);
if isempty(parse_and_order)
    parse_and_order=true(1,maxdepth);
end
for col=1:size(dimnames,2)
    nums=cell(size(dimnames{col},1),1);
    letters=cell(size(dimnames{col},1),1);
    for tok = 1:size(dimnames{col},1)
        nums{tok} = str2double(regexp(dimnames{col}{tok}, '\d+', 'match'));
        letters{tok} = regexp(dimnames{col}{tok}, '[A-z]+', 'match');
        if size(letters{tok,:},2)>1
            letters{tok}=letters{tok}(1);
        end
    end
   
      

    %Check if all strings have the same length
    if all(cellfun(@(x) isequal(size(letters{1}),size(x)),letters))
        dimchars=unique(vertcat(letters{:}),'stable');
    end
    if all(cellfun(@(x) isequal(size(nums{1}),size(x)),nums))
        nums_array= vertcat(nums{:});
    end
    
    try
    numerated_charspec(col) =  numel(dimchars)==numel(letters{1}) & numel(nums_array)>1 & parse_and_order(col);
    catch
        keyboard
    end
    dimchars=dimchars{1};
    nums_array=nums_array(:,1);
    if numerated_charspec(col)
        dimspecs{col,1}=unique(dimchars,'rows','stable');
        dimspecs{col,2}=unique(nums_array,'rows','stable');
%         strcat(dimspecs{4,1}{:})
    else
        dimspecs{col,1}=unique(dimnames{col},'stable');
        dimspecs{col,2}=[];
    end
    
end
pathListIdces=cell(size(pathlist,1),2);
treedepth=size(parsepath,2);
for row=1:size(parsepath,1)
    pathListIdces{row,1}=pathlist{row};
    idces=NaN(1,treedepth);
    for col=1:treedepth
        if size(dimspecs{col,1},1)==1
            idces(col)=str2double(strrep(regexp(parsepath{row,col}, [dimspecs{col,1} '\d+'], 'match'), dimspecs{col,1},''));
        else
            idces(col)=find(strcmp(dimspecs{col,1},parsepath{row,col}));
        end
    end
    pathListIdces{row,2}=idces;
end
[~, sorted_idces]=sortrows(vertcat(pathListIdces{:,2}));
pathListIdces=pathListIdces(sorted_idces,:);


trajectory_location=idSocial_buildTrajectoryLocationTreeFromList(pathListIdces);

%%
% trajectory_location = ...
%     idSocial_recursiveDelSubdivision(trajectory_location,superflous_subdivision);

