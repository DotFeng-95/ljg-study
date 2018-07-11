function cmp_mesh = get_cmp_mesh(initmesh, radius)
%%  generate valid meshgrid for computation

cmp_mesh.surface = initmesh.surface;
   cmp_mesh.elem = initmesh.elem;
   cmp_mesh.node = radius * initmesh.node;
   cmp_mesh.radius = radius;

%%%% computing surface area
%%%% Only for spherical surface 
%%%area = spherearea(cmp_mesh.node, cmp_mesh.elem);
%%%area_s = sum(area(:))

%%%% computing simplex area
cmp_mesh.area = simplexvolume(cmp_mesh.node, cmp_mesh.elem);
tol_area = sum(cmp_mesh.area(:));
cmp_mesh.tol_area = tol_area;

end
