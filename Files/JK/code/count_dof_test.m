load('../data/sph-cly-test1-5.mat');
rho = scft.rho;
node = init_mesh.node;
elem = init_mesh.elem;
n= countdots(node, elem, rho);
n 