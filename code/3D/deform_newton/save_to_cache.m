function [] = save_to_cache(kernel, params, results)
    %SAVE_TO_CACHE Save results to a dedicated cache file
    %
    % params: The struct containing the experiment parameters

    cache_file = sprintf('run_cache_%s_%s.mat', kernel, DataHash(params));
    save(cache_file, 'results', 'params');
end
