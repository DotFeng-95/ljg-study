function torussurface_test(type) 
unix('rm ./results/log.txt');
unix('rm ./results/iteration.dat');
unix('touch ./results/iteration.dat');

diary('./results/log.txt');
diary on;

sys_pmt.Nspecies = 2;
sys_pmt.Nblend   = 1;
sys_pmt.Nblock   = 2;
sys_pmt.Ndeg  	 = 100;

type = 4;
 
if type == 1
    %%%%  cylinder phase
    sys_pmt.fA  = 0.2;
    sys_pmt.chiAB = 0.25;
    R = 8;r = 1;nu = 320;nv = 40;
    initType = 1;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
elseif type == 2
    sys_pmt.fA  = 0.2;
    sys_pmt.chiAB = 0.25;
    R = 4;r = 2;nu = 320;nv = 40;
    initType = 1; 
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
elseif type == 3
    sys_pmt.fA  = 0.2;
    sys_pmt.chiAB = 0.25;
    R = 4;r = 2;nu = 320;nv = 40;
    initType = 2; 
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
elseif type == 4
    sys_pmt.fA  = 0.2;
    sys_pmt.chiAB = 0.25;
    R = 4/1.2;r = 2*1.2;nu = 320;nv = 40;
    initType = 1;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler'; 
    initrefine = 0;
    maxrefine = 1;
elseif type == 5
    %%%%  lamella phase
    sys_pmt.fA  = 0.5;
    sys_pmt.chiAB = 0.16;
    R = 8;r = 5;
    initType = 3 ;
    ntheta = 0;nphi = 8;nu = 160;nv = 100;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
elseif type == 6
    sys_pmt.fA  = 0.5;
    sys_pmt.chiAB = 0.16;
    R = 8;r = 5;
    initType = 3;
    ntheta = 5;nphi = 5;nu = 160;nv = 100;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
elseif type == 7
    sys_pmt.fA  = 0.5;
    sys_pmt.chiAB = 0.16;
    R = 8;r = 5;
    initType = 3; 
    ntheta = 6;nphi = 6;nu = 160;nv = 100;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
elseif type == 8
    sys_pmt.fA  = 0.5;
    sys_pmt.chiAB = 0.16;
    R = 8;r = 5;
    initType = 3; 
    ntheta = 7;nphi = 7;nu = 160;nv = 100;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
else type == 9
    sys_pmt.fA  = 0.5;
    sys_pmt.chiAB = 0.16;
    R = 8;r = 5;
    initType = 3; 
    ntheta = 8;nphi = 8;nu = 160;nv = 100;
    Radius = 1;
    Mthd.PDE = 'ImplicitCN';
    Mthd.Iter = 'Euler';
    initrefine = 0;
    maxrefine = 1;
end

chiN = sys_pmt.chiAB*sys_pmt.Ndeg;

surface = torussurface(R, r, nu, nv);

[node, elem, uv] = surface.initmesh();
ndof = size(node, 1);
mu = zeros(2, ndof);
w = zeros(2, ndof);

if initType == 1
    mu(1, :) = chiN*(-1 + 2*rand(1, ndof));
    mu(2, :) = chiN*(-1 + 2*rand(1, ndof));
elseif initType == 2
    mu(2, :) = chiN/2*(cos(3*uv(:, 1)).*cos(3*uv(:, 2)) + 0.5*cos(6*uv(:, 2)));
elseif initType == 3
    mu(2, :) = chiN*sin(ntheta*uv(:, 1) + nphi*uv(:, 2)); 
end


w(1,:) = mu(1,:) - mu(2,:);
w(2,:) = mu(1,:) + mu(2,:); 

init_mesh.surface = surface;
init_mesh.node = node; 
init_mesh.elem = elem;
%%%%% generate computational meshgrid  
cmp_mesh = get_cmp_mesh(init_mesh, Radius);
cmp_pmt = get_cmp_pmt(sys_pmt, init_mesh);
[scft, scftAux] = get_scft_data(sys_pmt, cmp_pmt); 

scft.w = w;
scft.mu(1,:) = 0.5*(scft.w(1,:)+scft.w(2,:));   
scft.mu(2,:) = 0.5*(scft.w(2,:)-scft.w(1,:)); 
scft.rho(1,:) = 0.5 + scft.mu(2,:)/chiN;
scft.rho(2,:) = 1.0 - scft.rho(1,:);

[cmp_mesh, cmp_pmt, scft] = diblockmain(init_mesh, cmp_mesh, cmp_pmt, scft, scftAux, sys_pmt, Mthd);

save(['torussurface', int2str(type), '.mat']) 

diary off;

