function plot_mode_out=idSocial_auxiliaries_setPlotMode(plot_mode,plot_mode_def)

if ~isempty(plot_mode_def)
    plot_mode_def_flds=fieldnames(plot_mode_def);
    no_fields_def=size(plot_mode_def_flds,1);
    plot_mode_out=plot_mode;
    
    for def_fld=1:no_fields_def
        if isfield(plot_mode_out,plot_mode_def_flds{def_fld})
            plot_mode_out.(plot_mode_def_flds{def_fld})=plot_mode_out.(plot_mode_def_flds{def_fld});
        elseif ~isfield(plot_mode_out,plot_mode_def_flds{def_fld}) && iscell(plot_mode_def.(plot_mode_def_flds{def_fld})) && ...
                numel(plot_mode_def.(plot_mode_def_flds{def_fld})) > 1 % Various options
            plot_mode_out.(plot_mode_def_flds{def_fld})=plot_mode_def.(plot_mode_def_flds{def_fld})(1);
        else
            plot_mode_out.(plot_mode_def_flds{def_fld})=plot_mode_def.(plot_mode_def_flds{def_fld});
        end
    end
else
    plot_mode_out = plot_mode;
end


