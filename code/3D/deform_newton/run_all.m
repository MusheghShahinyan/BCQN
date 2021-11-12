function [] = run_all()
%RUN_ALL Summary of this function goes here
%   Detailed explanation goes here

kernels = ["elephant", "shear_bar", "botijo", "armadillo"];

for i = 1:length(kernels)
    test(kernels(i))
end

end

