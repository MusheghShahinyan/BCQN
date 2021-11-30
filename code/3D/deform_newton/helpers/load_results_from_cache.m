function [success, results] = load_results_from_cache(kernel, re_run, params, params_index, save_to_workspace)
    %LOAD_RESULTS_FROM_CACHE Try to load results from a dedicated 
    % cache file 
    %
    % kernel: The name of the kernel 
    % re_run: The indices of experiments to re_run, as opposed to loading
    % params: The struct containing the experiment parameters
    % params_i: The param_group index

    global workspace_cache
    
    if nargin < 5
        save_to_workspace = 1;
    end

    if not(exist('workspace_cache','var')) 
        workspace_cache = struct();
    end
    
    cache_name = sprintf('run_cache_%s_%s', kernel, DataHash(params));
    cache_file = sprintf('cache_files/%s.mat', cache_name);
    

    if (isvector(re_run) && any(re_run == params_index)) || ...
            (isstring(re_run) && re_run == "all")
        results = [];
        success = false;
    elseif isfield(workspace_cache, cache_name)
        disp(['Using workspace cache for param_group(', num2str(params_index), ')']);
        results = workspace_cache.(cache_name);
        success = true;
    elseif isfile(cache_file)
        disp(['Using cache for param_group(', num2str(params_index), ')']);
        load(cache_file, 'results');
        
        if save_to_workspace
            workspace_cache.(cache_name) = results;
        end
        success = true;
    else
        results = [];
        success = false;
    end
end