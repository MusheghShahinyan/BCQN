function [success, results] = load_results_from_cache(kernel, re_run, params, params_index)
    %LOAD_RESULTS_FROM_CACHE Try to load results from a dedicated 
    % cache file 
    %
    % kernel: The name of the kernel 
    % re_run: The indices of experiments to re_run, as opposed to loading
    % params: The struct containing the experiment parameters
    % params_i: The param_group index

    global workspace_cache
    
    if not(exist('workspace_cache','var')) 
        workspace_cache = struct();
    end
    
    cache_name = sprintf('run_cache_%s_%s', kernel, DataHash(params));
    cache_file = sprintf('cache_files/%s.mat', cache_name);
    
    if isfield(workspace_cache, cache_name)
        disp(['Using workspace cache for param_group(', num2str(params_index), ')']);
        results = workspace_cache.(cache_name);
        success = true;
    elseif isfile(cache_file) && ...
       not( ...
            (isvector(re_run) && any(re_run == params_index)) || ...
            (isstring(re_run) && re_run == "all") ...
       )
        disp(['Using cache for param_group(', num2str(params_index), ')']);
        load(cache_file, 'results');
        
        workspace_cache.(cache_name) = results;
        success = true;
    else
        results = [];
        success = false;
    end
end