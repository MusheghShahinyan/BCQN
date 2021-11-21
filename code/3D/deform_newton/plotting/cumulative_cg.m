
function [] = cumulative_cg(param_group_results, param_groups, kernel)
%PLOT_ENERGY Summary of this function goes here
%   Detailed explanation goes here
figure; energy_axes = axes; hold(energy_axes, "on");

title(energy_axes, ['Energy Plot Cumulative CG', kernel,]);
xlabel(energy_axes, 'Cumulative CG Iterations');
ylabel(energy_axes, 'Energy Value');


for j = 1:length(param_groups)
    results = param_group_results{j};
    if isfield(results, 'num_iter')
        disp(param_groups(j).name);
        
        cumulative_iters = cumsum([0 results.num_iter']);

        stairs(energy_axes, cumulative_iters, results.energies, 'DisplayName', param_groups(j).name);
    end
    
    set(findall(energy_axes,'YAxisLocation','left'),'Yscale','log');

end

legend(energy_axes);
hold(energy_axes, "off");
saveas(energy_axes, strcat('figures/', kernel, '_energy.fig'));

end