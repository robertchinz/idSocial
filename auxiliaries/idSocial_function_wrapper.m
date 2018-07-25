function input_data=idSocial_function_wrapper(input_data, options, def_options_function, plot_mode, functionInfo, external_func_name)
if nargin < 6 || isempty(external_func_name)
    external_func_name = '';
end

if numel(input_data) == 1
    input_data=idSocial_function_wrapper_FlatInput(input_data, options, def_options_function, plot_mode, functionInfo, external_func_name);
else
    input_data=idSocial_function_wrapper_StructMultiDim(input_data, options, def_options_function, plot_mode, functionInfo, external_func_name);
end