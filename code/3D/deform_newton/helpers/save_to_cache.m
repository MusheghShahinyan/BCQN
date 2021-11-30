function [] = save_to_cache(kernel, params, results, save_to_workspace)
    %SAVE_TO_CACHE Save results to a dedicated cache file
    %
    % params: The struct containing the experiment parameters

    if nargin < 4
        save_to_workspace = 1;
    end

    global workspace_cache
    
    cache_name = sprintf('run_cache_%s_%s', kernel, DataHash(params));
    cache_file = sprintf('cache_files/%s.mat', cache_name);

    save(cache_file, 'results', 'params');
    if save_to_workspace
        workspace_cache.(cache_name) = results;
    end
end
