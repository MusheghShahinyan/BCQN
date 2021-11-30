% load bunny
global tet_mesh
load_input_mesh('/sdb/BCQN/code/3D/data/elephant');

loc = tet_mesh.v;
size(tet_mesh.t, 1)
faces = zeros(3, size(tet_mesh.t, 1)*4);

% Convert tets to triangles, transposes as well
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

%data = readply('bun_zipper.ply');
%loc = getfields(data.vertex,2,'x','y','z');
%faces = vertcat(data.face{:});

% make bunny upright (Z up)
%M = MRot3D([-90 0 0],1) * MScale3D(20);
%loc = loc * M(1:3,1:3);

f2 = vertcat(faces) - 1;
f2(1:16)

% mesh view, white
M = MRot3D([45 95 90],1) * MScale3D(20);
glViewer3D(loc * M(1:3,1:3), [121/255,139/255,179/255], vertcat(faces)' - 1);

% point cloud, auto color (scaled on Z)
% glViewer3D(loc);

% meshview, auto color
% glViewer3D(loc,[],faces);

% mesh view, confidence gray scale
% glViewer3D(loc,data.vertex.confidence,faces);



