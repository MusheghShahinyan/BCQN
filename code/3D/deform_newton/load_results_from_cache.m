function [success, results] = load_results_from_cache(kernel, re_run, params, params_index)
    %LOAD_RESULTS_FROM_CACHE Try to load results from a dedicated 
    % cache file 
    %
    % kernel: The name of the kernel 
    % re_run: The indices of experiments to re_run, as opposed to loading
    % params: The struct containing the experiment parameters
    % params_i: The param_group index

    cache_file = sprintf('run_cache_%s_%s.mat', kernel, DataHash(params));
    if isfile(cache_file) && ...
       not( ...
            (isvector(re_run) && any(re_run == params_index)) || ...
            (isstring(re_run) && re_run == "all") ...
       )
        disp(['Using cache for param_group(', num2str(params_index), ')']);
        load(cache_file, 'results');
        success = true;
    else
        results = [];
        success = false;
    end
end