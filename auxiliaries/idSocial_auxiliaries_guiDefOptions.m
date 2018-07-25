function [def_options, def_options_cell]=idSocial_auxiliaries_guiDefOptions(filename)

def_options=[];
def_options_cell=cell(1,5);
A=fileread(filename);
% expr='(?m)(?<=def_options.)(.*?)(?=(\r|%))';
expr='def_options\.[a-z_]+\=[^\r%]*(\r|%)\>';
B=regexp(A, expr,'match');


for k=1:size(B,2)
    if isempty(strfind(B{k},'act_method'))
        
        try
            eval([B{k}]);
        catch
            strsp=textscan(B{k},'%s','Delimiter','=');
            eval([strsp{1}{1} '=[]']);
        end
    end
end
if ~isempty(def_options)
    def_options_cell{:,1}= fieldnames(def_options); % Fieldnames
    def_options_cell{:,2}= struct2cell(def_options); % Struct values
    def_options_cell{:,3}= cellfun(@(x) class(x),def_options_cell{:,2},'UniformOutput',false); % Class of values
    def_options_cell{:,4}= cellfun(@(x) size(x),def_options_cell{:,2},'UniformOutput',false); % Size of values
    def_options_cell=cat(2,def_options_cell{:});
    
    numIdx = strcmp(def_options_cell(:,3),'double');
    charIdx = strcmp(def_options_cell(:,3),'char');
    cellIdx = strcmp(def_options_cell(:,3),'cell');
    structIdx = strcmp(def_options_cell(:,3),'struct');
    scalarIdx=cellfun(@(x) all(size(x)==1),def_options_cell(:,2));
    def_options_cell(numIdx & ~scalarIdx,3)=cellstr(repmat('numarray',sum(numIdx & ~scalarIdx,1),1));
  
    def_options_cell(charIdx,5)=def_options_cell(charIdx,2);
    def_options_cell(~charIdx & ~cellIdx & ~structIdx & scalarIdx,5)= cellfun(@(x) num2str(x), def_options_cell(~charIdx & ~cellIdx & ~structIdx & scalarIdx,2), 'UniformOutput',false  );
    
    for k=find(~scalarIdx & numIdx)'
        num2strarray='';
        for row=1:size(def_options_cell{k,2},1)
            for col=1:size(def_options_cell{k,2},2)
                if ischar(def_options_cell{k,2}(row,col))
                    %                     keyboard
                    num2strarray=[num2strarray '''' def_options_cell{k,2}(row,col) ''''  ', '];
                else
                    num2strarray=[num2strarray num2str(def_options_cell{k,2}(row,col)) ', '];
                end
            end
            num2strarray=num2strarray(1:end-2);
            num2strarray=[num2strarray '; '];
        end
        num2strarray=num2strarray(1:end-2);
        num2strarray=['[' num2strarray ']'];
        def_options_cell{k,5}=num2strarray;
    end
    
%     keyboard
    for k=find(cellIdx)'
        cell2strarray='';
        for row=1:size(def_options_cell{k,2},1)
            for col=1:size(def_options_cell{k,2},2)
                if ischar(def_options_cell{k,2}{row,col})
                    %                     keyboard
                    cell2strarray=[cell2strarray '''' def_options_cell{k,2}{row,col} ''''  ', '];
                else
                    cell2strarray=[cell2strarray num2str(def_options_cell{k,2}{row,col}) ', '];
                end
            end
            cell2strarray=cell2strarray(1:end-2);
            cell2strarray=[cell2strarray '; '];
        end
        cell2strarray=cell2strarray(1:end-2);
        cell2strarray=['{' cell2strarray '}'];
        def_options_cell{k,5}=cell2strarray;
    end
    
    for k=find(structIdx)'
        strct=def_options_cell{k,2};
%         strchar=evalc('disp(def_options_cell{:,2}{k})');
        fldnm=fieldnames(strct);
        no_opts=size(def_options_cell,1);
        reference_options=def_options_cell(k,:);%cellfun(@(x) x{k},def_options_cell,'UniformOutput',false);
        def_options_cell(k,:)=[];
     
        for fact=1:size(fldnm,1)
            struct2strarray='';
            
            newidx=no_opts+fact-1;
            def_options_cell{newidx,1}=[reference_options{1} '.' fldnm{fact}];
%             struct2strarray=[struct2strarray fldnm{fact} '='];
            act_fld=strct.(fldnm{fact});
            if  iscell(act_fld)
                def_options_cell{newidx,3}='struct';
                def_options_cell{newidx,2}=act_fld;
                def_options_cell{newidx,4}=size(act_fld);
                cell2strarray='';
                for row=1:size(act_fld,1)
                    for col=1:size(act_fld,2)
                        
                        if ischar(act_fld{row,col})
                            %                     keyboard
                            cell2strarray=[cell2strarray '''' act_fld{row,col} ''''  ', '];
                        elseif isnumeric(act_fld{row,col}) && all(size(act_fld{row,col})<2)
                            cell2strarray=[cell2strarray num2str(act_fld{row,col}) ', '];
                        elseif isnumeric(act_fld{row,col}) && any(size(act_fld{row,col})>=2)
                            num2strarray='';
                            for row2=1:size(act_fld{row,col},1)
                                for col2=1:size(act_fld{row,col},2)
                                    if ischar(act_fld{row,col}(row2,col2))
                                        %                     keyboard
                                        num2strarray=[num2strarray '''' act_fld{row,col}(row2,col2) ''''  ', '];
                                    else
                                        num2strarray=[num2strarray num2str(act_fld{row,col}(row2,col2)) ', '];
                                    end
                                end
                                num2strarray=num2strarray(1:end-2);
                                num2strarray=[num2strarray '; '];
                            end
                            num2strarray=num2strarray(1:end-2);
                            num2strarray=['[' num2strarray ']'];
                            cell2strarray=[cell2strarray num2strarray ', '];
                                                    
                        end
                    end
                    cell2strarray=cell2strarray(1:end-2);
                    cell2strarray=[cell2strarray '; '];
                    
                end
                cell2strarray=cell2strarray(1:end-2);
                cell2strarray=['{' cell2strarray '}'];
                
            end
            struct2strarray=[struct2strarray cell2strarray '; '];
            def_options_cell{newidx,5}=struct2strarray;
        end
        
    end
    
end
