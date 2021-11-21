function [] = plot_energy_grad(param_group_results, param_groups, kernel)
%PLOT_ENERGY_GRAD Summary of this function goes here
%   Detailed explanation goes here

figure; energy_axes_grad = axes; hold(energy_axes_grad, "on");
title(energy_axes_grad, ['Energy Plot Grad', kernel]);


for j = 1:length(param_groups)
    results = param_group_results{j};

    disp(param_groups(j).name);
    
    % Shift x axis to align with newton steps + starting point
    x_axis_energies = linspace(-1, length(results.energies) - 2, length(results.energies));
    % Shift x axis to align with newton steps
    x_axis = linspace(0, length(results.etas) - 1, length(results.etas));
    
    
    plot(energy_axes_grad, x_axis_energies, results.energies, 'DisplayName', param_groups(j).name);
    yyaxis(energy_axes_grad, 'right')
    plot(energy_axes_grad, x_axis, vecnorm(results.bs), 'DisplayName', param_groups(j).name);
    set(energy_axes_grad,'Yscale','log');
    yyaxis(energy_axes_grad, 'left')

end

xlabel(energy_axes_grad, 'Newton Iteration');
ylabel(energy_axes_grad, 'Energy Value');
set(energy_axes_grad,'Yscale','log');

legend(energy_axes_grad);

hold(energy_axes_grad, "off");

end

