function [ ] = opengl_render_result( u_x, frame )
%PLOT_RESULT Summary of this function goes here
%   Detailed explanation goes here

global tet_mesh K q ver_num

q_x = K * u_x + q;
r_q_x = reshape(q_x, 3, ver_num)';

faces = zeros(3, size(tet_mesh.t, 1)*4);

% Convert tets to triangles,
%  transposes at the same time in terms

for i = 1:size(tet_mesh.t, 1)
    triagnles_i = (i-1)*4 + 1;
    
    % Triangle 1
    faces(1, triagnles_i + 0) = tet_mesh.t(i, 1);
    faces(2, triagnles_i + 0) = tet_mesh.t(i, 2);
    faces(3, triagnles_i + 0) = tet_mesh.t(i, 3);
    
    % Triangle 2
    faces(1, triagnles_i + 1) = tet_mesh.t(i, 1);
    faces(2, triagnles_i + 1) = tet_mesh.t(i, 2);
    faces(3, triagnles_i + 1) = tet_mesh.t(i, 4);
    
    % Triangle 3
    faces(1, triagnles_i + 2) = tet_mesh.t(i, 1);
    faces(2, triagnles_i + 2) = tet_mesh.t(i, 3);
    faces(3, triagnles_i + 2) = tet_mesh.t(i, 4);
    
    % Triangle 4
    faces(1, triagnles_i + 3) = tet_mesh.t(i, 2);
    faces(2, triagnles_i + 3) = tet_mesh.t(i, 3);
    faces(3, triagnles_i + 3) = tet_mesh.t(i, 4);
end

M = MRot3D([45 95 90],1) * MScale3D(20);
% Convert faces to 0 indexing for opengl
glViewer3D(r_q_x * M(1:3,1:3), [121/255,139/255,179/255], vertcat(faces)' - 1);
end

