function [] = plot_energy(param_group_results, param_groups, kernel)
%PLOT_ENERGY Summary of this function goes here
%   Detailed explanation goes here
figure; energy_axes = axes; hold(energy_axes, "on");
grid(energy_axes, "on"); 

title(energy_axes, ['Energy Plot ', kernel,]);
xlabel(energy_axes, 'Newton Iteration');
ylabel(energy_axes, 'Energy Value');

disp("Plotting Energy => ");

colors = lines(length(param_groups));
energyOpts.LineWidth = 2;
iterOpts.LineWidth = 1;
iterOpts.LineStyle = '--';
iterOpts.Marker = 'o';


for j = 1:length(param_groups)
    results = param_group_results{j};
    disp(strcat("    ", param_groups(j).name));
    
    % Shift x axis to align with newton steps + starting point
    x_axis_energies = linspace(-1, length(results.energies) - 2, length(results.energies));
    % Shift x axis to align with newton steps
    x_axis = linspace(0, length(results.energies) - 2, length(results.energies) - 1);
    
    plot(energy_axes, x_axis_energies, results.energies, 'DisplayName', param_groups(j).name, 'Color', colors(j, :), energyOpts);

    if isfield(results, 'num_iter')
        yyaxis(energy_axes, 'right')
        ylim([0 250]);
        plot(energy_axes, x_axis, results.num_iter, 'DisplayName', strcat(param_groups(j).name, ' - iters'),  'Color', colors(j, :), 'MarkerFaceColor', colors(j, :), iterOpts);
        yyaxis(energy_axes, 'left');
    end
    
    set(findall(energy_axes,'YAxisLocation','left'),'Yscale','log');

end

% Only integer Iterations
energy_axes.XTick = unique(round(energy_axes.XTick));

legend(energy_axes);
hold(energy_axes, "off");
saveas(energy_axes, strcat('figures/', kernel, '_energy.fig'));

end

