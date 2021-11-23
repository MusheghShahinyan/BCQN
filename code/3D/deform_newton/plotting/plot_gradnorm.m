function [] = plot_gradnorm(param_group_results, param_groups, kernel)
%PLOT_GRADNORM Summary of this function goes here
%   Detailed explanation goes here
figure; gradnorm_axes = axes; hold(gradnorm_axes, "on");

title(gradnorm_axes, ['GradNorm Plot ', kernel,]);
xlabel(gradnorm_axes, 'Newton Iteration');
ylabel(gradnorm_axes, 'GradNorm');


disp("Plotting GradNorm => ");
max_gradnorm = 0;
for j = 1:length(param_groups)
    results = param_group_results{j};
    disp(strcat("    ", param_groups(j).name));
    
    % Shift x axis to align with newton steps + starting point
    % Shift x axis to align with newton steps
    x_axis = linspace(0, length(results.energies) - 2, length(results.energies) - 1);
    
    plot(gradnorm_axes, results.gradnorm, 'DisplayName', param_groups(j).name);
    max_gradnorm = max(max_gradnorm, max(results.gradnorm));
end

legend(gradnorm_axes);
ylim(gradnorm_axes, [0, max_gradnorm]);
set(gradnorm_axes, 'Yscale', 'log');

hold(gradnorm_axes, "off");
saveas(gradnorm_axes, strcat('figures/', kernel, '_gradnorm.fig'));
end

