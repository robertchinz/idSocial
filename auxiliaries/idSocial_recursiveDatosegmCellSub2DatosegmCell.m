function [patterncell datosegmcell] = idSocial_recursiveDatosegmCellSub2DatosegmCell(patterncell,datosegmcell)
% This function takes the cell 'patterncell', which is a 'cell of cell' of the form patterncell{1}=cell(..,..);
% patterncell{1}{1}=cell(...,...); etc. and copies information
% from 'datosegmcell', which has the samme structure but can
% have smaller depth downstream in a way that if
% datosegmcell{1}{1} = options, patterncell{1}{1}{:}=options.
% ATTENTION: The first output patterncell returns the new
% datosegmcell!

if nargin<2
    datosegmcell=[];
end

if isa(patterncell,'cell')%all(cellfun(@(x) isa(x,'cell'),patterncell))
    for k = 1:numel(patterncell)
        if numel(datosegmcell)>=k && isa(datosegmcell,'cell')
            
            [patterncell{k} datosegmcell{k}] = idSocial_recursiveDatosegmCellSub2DatosegmCell(patterncell{k},datosegmcell{k});
            
        else
            
            [patterncell{k} datosegmcell] = idSocial_recursiveDatosegmCellSub2DatosegmCell(patterncell{k},[]);
            
        end
        
        
        
    end
else
    
    
    patterncell=datosegmcell;
    return
end
end


