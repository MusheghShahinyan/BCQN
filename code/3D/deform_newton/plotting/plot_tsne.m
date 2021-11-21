function [] = plot_tsne(param_group_results, param_groups, kernel, u_n)
%PLOT_TSNE Summary of this function goes here
%   Detailed explanation goes here

figure; tsne_axes = axes; hold(tsne_axes, "on");

total_guesses = 0;
for j = 1:length(param_groups)
    results = param_group_results{j};
    total_guesses = total_guesses + size(results.guesses, 1);
end

tsne_X = zeros(total_guesses, length(u_n));
cur_tsne_count = 1; 

for j = 1:length(param_groups)        
    results = param_group_results{j};
    tsne_X(cur_tsne_count:(cur_tsne_count+size(results.guesses, 1)-1), :) = results.guesses;
    cur_tsne_count = cur_tsne_count + size(results.guesses, 1);
end


tsne_Y = tsne(tsne_X);
disp("Plotting TSNE");
gscatter(tsne_axes, tsne_Y(:, 1), tsne_Y(:, 2));
legend(tsne_axes);
title(tsne_axes, "Search Path Embeddign");
hold(tsne_axes, "off");

end

