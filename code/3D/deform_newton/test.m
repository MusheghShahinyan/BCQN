function [] = test(kernel)
addpath(genpath('../'))

if kernel == "elephant"
    example_elephant_init 
elseif kernel == "shear_bar"
    example_shear_bar_init
elseif kernel == "botijo"
    example_botijo_init
elseif kernel == "armadillo"
    example_armadillo_init
end

% example_elephant_init
% example_large_init
% example_shear_bar_init
% example_twist_bar_init
% example_stretch_bar_init
% example_rand_mesh_init
% example_zero_mesh_init
% example_armadillo_init
% example_dancer_init
% example_botijo_init
% example_dilo_test
% example_homer_init
% example_horse_init
% example_cube_init
% example_wrench_init
% example_statue_init
% example_armadillo_dance_init
% example_homer_bend_init


preconditioner = "incomplete_LU";
% preconditioner = 'diagonal'; 

pcg_parameters = struct( ...
    'tol', 1e-4, ...
    'restart', 100, ...
    'maxit', 100, ...
    'use_warm_start', false, ...
    'energy_tol', 1e-4);

use_direct = false; 

newton_solver(u_n, preconditioner, pcg_parameters, use_direct, kernel);

end

