function [] = plot_energy_grad(param_group_results, param_groups, kernel)
%PLOT_ENERGY_GRAD Summary of this function goes here
%   Detailed explanation goes here

figure; energy_axes_grad = axes; hold(energy_axes_grad, "on");
grid(energy_axes_grad, "on");
title(energy_axes_grad, ['Energy Gradient Plot', kernel]);

disp("Plotting Ennergy Gradients => ");

colors = lines(length(param_groups));
energyOpts.LineWidth = 2;

for j = 1:length(param_groups)
    results = param_group_results{j};
    disp(strcat("    ", param_groups(j).name));
    
    x_axis = linspace(0, length(results.gradnorms) - 1, length(results.gradnorms));
    plot(energy_axes_grad, x_axis, results.gradnorms, 'DisplayName', param_groups(j).name, 'Color', colors(j, :), energyOpts);
end

set(energy_axes_grad,'Yscale','log');

% Only integer Iterations
energy_axes_grad.XTick = unique(round(energy_axes_grad.XTick));

xlabel(energy_axes_grad, 'Newton Iteration');
ylabel(energy_axes_grad, 'Energy Value');
set(energy_axes_grad,'Yscale','log');

legend(energy_axes_grad);

hold(energy_axes_grad, "off");

end

