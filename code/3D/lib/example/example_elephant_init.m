
global ver_num tet_mesh dirichlet en_type cmin cmax K q

% colormap(jet);

en_type = 'arap';

load_input_mesh('../data/elephant');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = reshape(tet_mesh.v(:, :)', ver_num * 3, 1);

% q = rand(ver_num * 3, 1);

anchor1 = find(tet_mesh.v(:, 2) > 0.45);
anchor2 = find(tet_mesh.v(:, 3) < -0.2);
anchor3 = find(tet_mesh.v(:, 3) < -0.4);

mask = [anchor1; anchor2; anchor3];

q(3 * anchor1 - 0) = q(3 * anchor1 - 0) - 0.6;
q(3 * anchor1 - 1) = q(3 * anchor1 - 1) + 0.1;

dirichlet = zeros(3 * ver_num, 1);
dirichlet([3 * mask - 2, 3 * mask - 1, 3 * mask]) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


basis = find(dirichlet == 0);
num_basis = size(basis, 1);
u_n = zeros(num_basis, 1);
k_val = ones(num_basis, 1); 

K = sparse(basis, [1:num_basis]', k_val, 3 * ver_num, num_basis);


global_cache_claim


cmax = 1;


if strcmp(en_type, 'arap')
   
    cmin = 0;
    
end

if strcmp(en_type, 'mips')
   
    cmin = 2;
    
end

if strcmp(en_type, 'iso')
   
    cmin = 4;
    
end

if strcmp(en_type, 'amips')
   
    cmin = exp(amips_s * 2);
    
end


if strcmp(en_type, 'conf')
   
    cmin = 1;
    
end

