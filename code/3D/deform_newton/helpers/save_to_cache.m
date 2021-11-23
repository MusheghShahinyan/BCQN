function [] = save_to_cache(kernel, params, results)
    %SAVE_TO_CACHE Save results to a dedicated cache file
    %
    % params: The struct containing the experiment parameters

    global workspace_cache
    
    cache_name = sprintf('run_cache_%s_%s', kernel, DataHash(params));
    cache_file = sprintf('cache_files/%s.mat', cache_name);

    save(cache_file, 'results', 'params');
    workspace_cache.(cache_name) = results;
end
