function [success, u_n] = initialize_kernel(kernel)
    %INITIALIZE_KERNEL Set up the experiment by initializing global 
    % variables
    
    success = false;
    u_n = [];
    switch kernel 
        case "elephant"
            example_elephant_init; success = true;
        case "shear_bar"
            example_shear_bar_init; success = true;
        case "botijo"
            example_botijo_init; success = true;
        case "armadillo"
            example_armadillo_init; success = true;
        case "wrench"
            example_wrench_init; success = true;
        case "homer" 
            example_homer_bend_init; success = true;
        case "horse"
            example_horse_init; success = true;
        case "cube"  
            example_cube_init; success = true;
        case "dancer"    
            example_dancer_init; success = true;
        case "twist"  
            example_twist_bar_init; success = true;
        case "santa"  
            example_santa_init; success = true;
        case "statue"  
            example_statue_init; success = true;
    end
    
    if not(success)
        disp(["Failed to initialize the kernel: ", kernel]);
    end
end

