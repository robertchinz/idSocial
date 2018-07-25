function [good_files, idxOut]= idSocial_createFileTreeFromPath(PathName,ids,countable,cfilter,sfilters) 
% % PathName = 'G:\Ontogeny\data';
if nargin<2 || isempty(ids)
    ids=[];
end

if nargin<3 || isempty(countable)
    countable=[];
end
if nargin<4 || isempty(cfilter)
    cfilter=[];
end
if nargin<5 || isempty(sfilters)
    sfilters=[];
end

% if ~isempty(ids)
%     ids(cellfun(@(x) isempty(x), ids))= [];
% end

filelist = idSocial_recursiveDir(PathName,'*.mat');

% ids = {'N2','day'};
% sfilters = {'WT','trajectories.mat'};
% countable = [];%{false,false,true,false};
% countable = {false,true,false};
% cfilter = {[],[5 25],[]};


if isempty(countable)
    countable = repmat({false},1,numel(ids));
end
if ~isempty(cfilter) && ~isempty(countable)
    for k = 1:numel(countable)
        if ~countable{k} && numel(cfilter)>=k
            cfilter{k} = [];
        end
    end
end

not_isempty_id =cellfun(@(x) ~isempty(x),ids);
% Regexps:
% one integer: .$  (e.g., 'day.$' for day9 but not day10)

good_file = false(size(filelist,1),1);
idIdx = NaN(size(filelist,1),numel(ids));
realIds = cell(size(filelist,1),numel(ids));
idisfile = false(size(filelist,1),1);
for k = 1:size(filelist,1)
    % First check simple filter
    simpleFilter = true;
    if ~isempty(sfilters)
        simpleFilter = sum(~cellfun(@isempty,regexp(filelist{k},sfilters)))==numel(sfilters);
    end
    if simpleFilter
        fileprts = textscan(filelist{k},'%s','Delimiter',['\' filesep]);
        fileprts = fileprts{1};
        
        if numel(fileprts)-1>=numel(ids) % There must be enough subfolders to satisfy identifiers; -1 for drive letter
            for k2 = 1:numel(ids)
                %             boolIdx = ~cellfun(@isempty,strfind(fileprts,ids{k2}));
                boolIdx = ~cellfun(@isempty,regexp(fileprts,ids{k2}));
                if sum(boolIdx)==1
                    idIdx(k,k2) = find(boolIdx);
                    realIds(k,k2) = fileprts(boolIdx);
                    if idIdx(k,k2)==numel(boolIdx)
                        idisfile(k)=true;
                    end
                end
            end
            
        end
        
        good_file(k) = ~any(isnan(idIdx(k,not_isempty_id)));
    end
end
%% Find unique identifiers

realIds2 = realIds(good_file,:);
idisfile2 = idisfile(good_file);
uniqueRealIds = cell(numel(ids),1);
for k2 = 1:numel(ids)
    uniqArray =realIds2(:,k2);
    if ~all(cellfun(@(x) isempty(x),uniqArray))
        uniqueRealIds{k2} = unique(uniqArray);
    end
end

% Check if one of the ids correspond to file, not folder:
idisfile_all = all(idisfile2);

% See if we can put it in <= 3 subgroups
% Check if subgroups would be <=3 (that is all idSocial can handle right
% now
no_uniqueIds = cellfun(@(x) numel(x), uniqueRealIds);
if sum(no_uniqueIds>1)>3 || idisfile_all && sum(no_uniqueIds>1)==4
    errordlg('Sorry, idSocial can only handle 2 levels at the moment.')
end

good_ids = true(1,numel(no_uniqueIds));
% Clean unneccessary partitions

good_ids(no_uniqueIds==1)=false;
cleanUniqueRealIds = uniqueRealIds(good_ids);
countable = countable(good_ids);
cfilter = cfilter(good_ids);
ids = ids(good_ids);


%% Sort files into subgroups, etc.

good_files = filelist(good_file);
idxOut = ones(size(good_files,1),3);
for idCount = numel(cleanUniqueRealIds):-1:1
    for id1=1:numel(cleanUniqueRealIds{idCount})
        thisIdx = cellfun(@(x) ~isempty(x),strfind(good_files,cleanUniqueRealIds{idCount}{id1}));
        idsAct=id1;
        % If set as countable:
        if countable{idCount}
            try
                idsAct = str2double(strrep(cleanUniqueRealIds{idCount}{id1},ids{idCount},''));
                if ~isnan(idsAct) && isfinite(idsAct) && isnumeric(idsAct)
                    if ~isempty(cfilter) && ~isempty(cfilter{idCount}) && numel(cfilter{idCount})==2
                        if idsAct<cfilter{idCount}(1) || idsAct>cfilter{idCount}(2);
                            idsAct = -1;
                        end
                        
                    end
                else
                    idsAct=id1;
                    
                end
            catch
                warning([mfilename ': Sorry, could not extract subgroup count.']);
            end
        end
        idxOut(thisIdx,idCount) = idsAct;
    end
end

not_filtered = ~any(idxOut<0,2);
good_files = good_files(not_filtered);
idxOut = idxOut(not_filtered,:);


uniqueCombs = unique(idxOut(:,all(~isnan(idxOut),1)),'rows');

% Recount indices for non-countables in case there are index-gaps after
% filtering

for id = 1:size(ids,2)
    if ~countable{id}
        uniqueCombsUnique = unique(uniqueCombs(:,id));
        uniqueCombsCount = 1:numel(uniqueCombsUnique);
        for uid = 1:numel(uniqueCombsCount)
            idxOut(idxOut(:,id)==uniqueCombsUnique(uid),id) = uniqueCombsCount(uid);
        end
       
    end
end


for ids = 1:size(uniqueCombs,1)
    membsIdx = ismember(idxOut(:,all(~isnan(idxOut),1)),uniqueCombs(ids,:),'rows');
    trialIdx = cumsum(membsIdx);
    trialIdx = trialIdx(membsIdx);
    idxOut(membsIdx,3) = trialIdx;
end


% for ids = numel(cleanUniqueRealIds):-1:1
%     for id1=1:numel(cleanUniqueRealIds{ids})
%         for k=1:size(good_files,1)
%             
%             if~isempty(strfind(good_files{k},cleanUniqueRealIds{ids}{id1}));
%                 idxOut(k,ids) = id1;
%             end
%         end
%     end
% end

% Check if ids are unique

all_array = [num2cell(idxOut) good_files];
all_array = sortrows(all_array,[1 2 3]);

idxOut = cell2mat(all_array(:,1:3));
good_files = all_array(:,4);
%%
for k=1:size(good_files,1)
    disp(sprintf('%s %d %d %d', good_files{k},idxOut(k,1),idxOut(k,2),idxOut(k,3)))
end


% gui.openFolderFileChooserList = filelist;
% guidata(fh,gui);