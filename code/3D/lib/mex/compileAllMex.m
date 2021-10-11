function compileAllMex() 
    eigenFolder = '/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3';
    iglFolder = '/Users/musheghshahinyan/Projects/libigl/include';
    tbbFolder = '/opt/homebrew/Cellar/tbb/2021.4.0/include';
    tbbLib = '/opt/homebrew/Cellar/tbb/2021.4.0/lib';
	% tbbFolder = 'C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.0.109\windows\tbb\include';
	% eigenFolder = 'C:\Library\eigen';
	% iglFolder = 'C:\Library\libigl\include';

	mex(['-I',eigenFolder],'energy_hessian_mex.cpp');
	mex(['-I',tbbFolder], ['-I',eigenFolder], ['-L', tbbLib], ['-l', 'tbb'],'lcp_solve_tbb_mex.cpp');
	% mex(['-I',tbbFolder], ['-I',eigenFolder],'lcp_solve_tbb_mex.cpp');
	mex(['-I',eigenFolder],'lcp_solve_eigen_mex.cpp');
	mex(['-I',eigenFolder],'load_tet_mex.cpp');
	mex(['-I',eigenFolder],'save_tet_mex.cpp');
	mex('-largeArrayDims',['-I',eigenFolder],'precompute_mex.cpp');
	mex('-largeArrayDims',['-I',eigenFolder],'get_feasible_steps_mex.cpp');
	mex(['-I',eigenFolder],'energy_value_mex.cpp');
	mex(['-I',eigenFolder],'energy_color_mex.cpp');
	mex(['-I',eigenFolder],'grad_function_mex.cpp');
end

