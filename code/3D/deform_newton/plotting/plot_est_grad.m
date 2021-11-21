function [] = plot_est_grad(param_group_results, param_groups, kernel)
%PLOT_EST_GRAD Summary of this function goes here
%   Detailed explanation goes here

figure; est_grad = axes; hold(est_grad, "on");
title(est_grad, ['Est. Grad', kernel]);

est_grads = zeros(length(param_groups), 100);
act_grads = zeros(length(param_groups), 100);
labels = cell(1, length(param_groups)*2);
max_timestep = 0;

for j = 1:length(param_groups)
    results = param_group_results{j};
    
    % size(results.grads_pre_ls)
    estimated_grad = zeros(size(results.grads_pre_ls));
    
    for i = 1:length(results.grads_pre_ls(1,:))
        H = results.hessians_pre_ls{i};
        grad = results.grads_pre_ls(:, i);
        p = results.search_directions(:, i);
        estimated_grad(:, i) = grad + H*p;
    end
    
    labels{j} = param_groups(j).name + " est. grad";
    labels{length(param_groups) + j} = ...
        param_groups(j).name + " act. grad";
    
    est_grads(j, 1:length(estimated_grad(1,:))) = ...
        vecnorm(estimated_grad);
    act_grads(j, 1:length(results.grads_pre_ls(1,:))) = ...
        vecnorm(results.grads_pre_ls);
    
    disp(param_groups(j).name);

    if max(max_timestep, length(results.bs(1,:)) > max_timestep)
        max_timestep = max(max_timestep, length(results.bs(1,:)));
    end

end

% plot grad differences
est_grads = est_grads(:, 1:max_timestep);
act_grads = act_grads(:, 1:max_timestep);

% shift est. grads since they estimate the gradient of the next
%   timestep
est_grads = [zeros(length(est_grads(:,1)), 1) est_grads];
act_grads = [act_grads zeros(length(est_grads(:,1)), 1)];

len = length(act_grads(1,:));
x_axis = linspace(0, len - 1, len);
bar(est_grad, x_axis, [est_grads; act_grads])
set(est_grad,'Yscale','log');
legend(est_grad, labels);

hold(est_grad, "off");
end

