function [] = run_all(make_plots)
%RUN_ALL Summary of this function goes here
%   Detailed explanation goes here

kernels = ["shear_bar", "elephant", "botijo", "armadillo"];
% kernels = ["shear_bar", "shear_bar"];

energies = cell(length(kernels), 1);
for i = 1:length(kernels)
    test(kernels(i), make_plots);
end

end

