function init_mesh = get_init_mesh(surf_kind)
%%  generate initial meshgrid

if strcmp(surf_kind, 'spheresurface')
%%%% spherical surface
	surface = spheresurface();
	[node,elem] = surface.initmesh();
	N = size(node, 1);
	%%%% mesh refinement
	for i = 1:4
		N = size(node, 1);
		[node, elem] = uniformrefine(node,elem);
		node(N+1:end,:) = surface.project(node(N+1:end, :));
	end
    ndof = size(node, 1);
    mu = zeros(2, ndof);
    w = zeros(2, ndof);
    mu(2, 1:12) = 1;
    w(1,:) = mu(1,:) - mu(2,:);
    w(2,:) = mu(1,:) + mu(2,:);
elseif strcmp(surf_kind, 'torus')
    R = 8;
    r = 1;
    chiN = 25;
    initType = 1;
    R = 8;
    r = 1;
    nu = 320;
    nv = 40;

    ntheta = 5;
    nphi = 5;

    surface = torussurface(R, r, nu, nv);
    [node, elem, uv] = surface.initmesh();
    ndof = size(node, 1);
    mu = zeros(2, ndof);
    w = zeros(2, ndof);
    if initType == 1
        mu(1, :) = chiN*(-1 + 2*rand(1, ndof));
        mu(2, :) = chiN*(-1 + 2*rand(1, ndof));
    elseif initType == 2
        mu(2, :) = chiN/2*(cos(3*uv(:, 1))*cos(3*uv(:, 2)) + 0.5*cos(6*uv(:, 2)));
    elseif initType == 3
        mu(2, :) = chiN*sin(ntheta*uv(:, 1) + nphi*uv(:, 2));
    end

    w(1,:) = mu(1,:) - mu(2,:);
    w(2,:) = mu(1,:) + mu(2,:);
elseif strcmp(surf_kind, 'orthocircle')
	surface = orthocirclesurface();
	[node,elem] = surface.initmesh('./meshdata/orthocirclesurfaceopt.mat');
	%[node,elem] = optsurfacemesh(node,elem,surface, 1.5, 200);
	%save('meshdata/orthocirclesurfaceopt.mat', 'node', 'elem')

elseif strcmp(surf_kind, 'heartsurface')
	%%%%%  some problems of the following lines
     surface = heartsurface();
     [node,elem] = surface.initmesh('./meshdata/heartsurfaceopt.mat');

	%surface = heartsurface();
	%[node,elem] = surface.initmesh();


	surface = heartsurface();
	[node,elem] = surface.initmesh();

elseif strcmp(surf_kind, 'ellipsoid')
	%%%%%  some problems of the following lines
	surface = ellipsoidsurface();
	[node,elem] = surface.initmesh('./meshdata/ellipsoid3394.mat');
elseif strcmp(surf_kind,'doubletorus')
	%%%%%  some problems of the following lines
	surface = doubletorus();
	[node,elem] = surface.initmesh('./meshdata/doubletorussurface.mat');	

end 

   init_mesh.node = node;
   init_mesh.elem = elem; 
   init_mesh.surface = surface;
   init_mesh.w = w;
end
