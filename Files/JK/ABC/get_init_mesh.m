function init_mesh = get_init_mesh(surf_kind)
%%  generate initial sphere meshgrid

if strcmp(surf_kind, 'spheresurface')
%%%% spherical surface
	surface = spheresurface();
	[node,elem] = surface.initmesh();
	%%%% mesh refinement
	for i = 1:4
		N = size(node, 1);
		[node, elem] = uniformrefine(node,elem);
		node(N+1:end,:) = surface.project(node(N+1:end, :));
    end
    
end 
   init_mesh.node = node;
   init_mesh.elem = elem; 
   init_mesh.surface = surface;
end
