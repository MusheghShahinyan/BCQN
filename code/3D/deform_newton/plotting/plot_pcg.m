function [] = plot_pcg(param_group_results, param_groups, kernel)
%PLOT_PCG Summary of this function goes here
%   Detailed explanation goes here
disp("Plotting PCG")

for j = 1:length(param_groups)
    results = param_group_results{j};
    pcg_parameters = param_groups(j).pcg_parameters;

    adaptive_pcg = not( ...
        param_groups(j).use_direct) && ...
        isfield(pcg_parameters, 'energy_tol') && ...
        isfield(pcg_parameters, 'line_check_jump' ...
    );
    
    if adaptive_pcg
 
        figure; pcg_1 = subplot(3,1,1); pcg_2 = subplot(3,1,2); pcg_3 = subplot(3,1,3); 
        hold(pcg_1, "on"); hold(pcg_2, "on"); hold(pcg_3, "on"); 

        for iter = 1:length(results.resvecs)
            energyvec = results.energyvecs{iter};
            plot(pcg_1, (energyvec - min(energyvec)) / (energyvec(1) - min(energyvec)), 'DisplayName', ['iter ', num2str(iter - 1)]);
            plot(pcg_2, results.resvecs{iter}, 'DisplayName', ['iter ', num2str(iter - 1), 'res']);
            plot(pcg_3, results.anglesvecs{iter}, 'DisplayName', ['iter ', num2str(iter - 1), 'res']);
        end
        
        set(pcg_1,'Yscale','linear');
        legend(pcg_1);
        title(pcg_1, "Energy");
        hold(pcg_1, "off");
        
        set(pcg_2,'Yscale','log');
        legend(pcg_2);
        title(pcg_2, "Residuals")
        hold(pcg_2, "off");
        
        set(pcg_3,'Yscale','linear');
        legend(pcg_3);
        title(pcg_3, "Angles (w.r.t -grad)")
        hold(pcg_3, "off");
    end

end

end

