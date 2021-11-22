function [] = plot_energy(param_group_results, param_groups, kernel)
%PLOT_ENERGY Summary of this function goes here
%   Detailed explanation goes here
figure; energy_axes = axes; hold(energy_axes, "on");
figure; gradnorm_axes = axes; hold(gradnorm_axes, "on");

title(energy_axes, ['Energy Plot ', kernel,]);
xlabel(energy_axes, 'Newton Iteration');
ylabel(energy_axes, 'Energy Value');

title(gradnorm_axes, ['GradNorm Plot ', kernel,]);
xlabel(gradnorm_axes, 'Newton Iteration');
ylabel(gradnorm_axes, 'GradNorm');


disp("Plotting Energy => ");
mmm = 0;
for j = 1:length(param_groups)
    results = param_group_results{j};
    disp(strcat("    ", param_groups(j).name));
    
    % Shift x axis to align with newton steps + starting point
    x_axis_energies = linspace(-1, length(results.energies) - 2, length(results.energies));
    % Shift x axis to align with newton steps
    x_axis = linspace(0, length(results.energies) - 2, length(results.energies) - 1);
    
    plot(energy_axes, x_axis_energies, results.energies, 'DisplayName', param_groups(j).name);
    plot(gradnorm_axes, results.gradnorm, 'DisplayName', param_groups(j).name);
    mmm = max(mmm, max(results.gradnorm));

    if isfield(results, 'num_iter')
        yyaxis(energy_axes, 'right')
        ylim([0 250]);
        plot(energy_axes, x_axis, results.num_iter, 'DisplayName', strcat(param_groups(j).name, '- iters'));
        yyaxis(energy_axes, 'left');
    end
    
    set(findall(energy_axes,'YAxisLocation','left'),'Yscale','log');

end

legend(energy_axes);
hold(energy_axes, "off");
saveas(energy_axes, strcat('figures/', kernel, '_energy.fig'));

legend(gradnorm_axes);
ylim(gradnorm_axes, [0, mmm]);
set(gradnorm_axes, 'Yscale', 'log');

hold(gradnorm_axes, "off");
saveas(gradnorm_axes, strcat('figures/', kernel, '_gradnorm.fig'));
end

