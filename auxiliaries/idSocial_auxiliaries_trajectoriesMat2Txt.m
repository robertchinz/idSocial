% function idSocial_auxiliaries_trajectoriesMat2Txt(tr_location)
addprobs = true;
tr = load(tr_location);
probtr = tr.probtrajectories;
tr = tr.trajectories;

no_fish = size(tr,2);
no_frames = size(tr,1);
no_dims = size(tr,3);
if addprobs
    no_probs = 1;
else
    no_probs = 0;
end

out = [tr_location(1:end-4) '.txt'];

fileID = fopen(out,'wt');

formatSpecH = repmat('%s\t',[1 no_fish*(no_dims+no_probs)]);
formatSpecH  = formatSpecH(1:end-2);
formatSpecH = [formatSpecH '\n'];
headerstring = cell(1,no_fish*(no_dims+no_probs));
colCount=1;
for ff = 1:no_fish
    for dm = 1:(no_dims+no_probs)
       switch dm
           case 1
               dimStr = 'X';
           case 2
               dimStr = 'Y';
           case 3
               dimStr = 'Z';
       end
       if dm == no_dims+no_probs && no_probs == 1
           headerstring{:,colCount} = ['P' num2str(ff)];
       else
            headerstring{:,colCount} = [dimStr num2str(ff)];
       end
        colCount = colCount + 1;
    end
end

fprintf(fileID,formatSpecH,headerstring{:});

formatSpec = repmat('%.2f\t',[1 no_fish*(no_dims+no_probs)]);
formatSpec  = formatSpec(1:end-2);
formatSpec = [formatSpec '\n'];

out_array = NaN(no_frames,no_fish*(no_dims+no_probs));
colCount=1;
for ff = 1:no_fish
    for dm = 1:no_dims+no_probs
        if dm == no_dims+no_probs && no_probs == 1
            out_array(:,colCount) = probtr(:,ff);
        else
            out_array(:,colCount) = tr(:,ff,dm);
        end
        colCount = colCount + 1;
    end
end

fprintf(fileID,formatSpec,out_array);
fclose(fileID);