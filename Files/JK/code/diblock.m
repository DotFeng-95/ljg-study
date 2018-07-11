clear all;
close all
clc;

format long;

unix('rm ./results/log.txt');
diary('./results/log.txt');
diary on;

addpath(genpath(pwd),'-begin');
savepath;

char saveName;

%%%% set model parameters 
sys_pmt = get_sys_pmt();

%%%% set initial surface meshgrid
surf_kind = 'spheresurface';
%surf_kind = 'orthocircle';
%surf_kind = 'heartsurface';
%surf_kind = 'ellipsoid';
%surf_kind = 'torus'

init_mesh = get_init_mesh(surf_kind);

cmp_pmt = get_cmp_pmt(sys_pmt, init_mesh);
[scft, scftAux] = get_scft_data(sys_pmt, cmp_pmt);

%%%%% set flag for numerical schemes
%%%%% the method to solve PDEs
%Mthd.PDE = 'RK4';  %%% unstable
%Mthd.PDE = 'ExplictEuler';
%Mthd.PDE = 'ImplicitEuler';
Mthd.PDE = 'CN';

%%%%% the method to update fields
Mthd.Iter = 'Euler';

%%%%  radius of surface
%%%%  lamella phase
%Radius = 8.0;
%%%%  cylinder phase oin spheresurface 
Radius = 3.65;

%%%%  cylinder phases on torus surface
%Radius = 1;
%%%%% generate computational meshgrid
cmp_mesh = get_cmp_mesh(init_mesh, Radius);

nspecies = sys_pmt.Nspecies;
	chiN = sys_pmt.Ndeg*sys_pmt.chiAB;
	nt    = cmp_pmt.Nt;
	ndof  = cmp_pmt.Ndof;
	nelem = cmp_pmt.Nelem;
	tol   = cmp_pmt.TOL;

fprintf('\n\n============== START RUNNING ... ==============\n\n\n');

lid = fopen('./results/iteration.dat', 'w');
fprintf(lid, '');
fclose(lid);

fprintf('========== Model Parameters ==========\n\n');
fprintf('chiN = %.3f, fA = %.3f\n\n', chiN, sys_pmt.fA);

fprintf('========== Discretization Points ==========\n\n');
fprintf('Dofs: [Nt, Nelem, Ndof] = [%d, %d, %d], PDEsolver: %s\n\n\n', nt, nelem, ndof, Mthd.PDE);

%%%  w: nspecies x ndof
choice_Field = 'Mu_Field';

%%% IMPORTANCE: obtain initial fields 
%scft.w = initial_field(scft.w, cmp_mesh, sys_pmt, cmp_pmt);
%%% IMPORTANCE: obtain initial fields 
scft.w = init_mesh.w;

scft.mu(1,:) = 0.5*(scft.w(1,:)+scft.w(2,:));
scft.mu(2,:) = 0.5*(scft.w(2,:)-scft.w(1,:));

%%%%  Initial values TEST
scft.rho(1,:) = 0.5 + scft.mu(2,:)/chiN;
scft.rho(2,:) = 1.0 - scft.rho(1,:);

figure(1)
showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));
saveName = sprintf('./results/phiA.%d', 0);
title(saveName);
set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
f = getframe(gcf);
imwrite(f.cdata,  [saveName, '.png']);
%pause();

%%%% only for explicit scheme of solving PDE
invMat = assembleMatrix(cmp_mesh, Mthd);
Mthd.invMat = invMat;
%%%% only for explicit scheme of solving PDE

dtol_energy = inf;
res = inf;
ITER = 0;

t1loop = clock;
[energy, scft] = evlSaddle(cmp_pmt, cmp_pmt.Maxiter, sys_pmt, cmp_mesh, scft, scftAux, Mthd);
old_energy = energy;

fprintf('\n &&&&&& ITER %d : Radius = %.15e, energy = %.15e\n\n', ITER, Radius, energy);

while res > cmp_pmt.TOL
	ITER = ITER + 1;

	Radius = adaptSurface(Radius, scft, cmp_pmt, init_mesh, Mthd, 1.0e-3);

	%%%%% generate computational meshgrid
	cmp_mesh = get_cmp_mesh(init_mesh, Radius);
	[energy, scft] = evlSaddle(cmp_pmt, cmp_pmt.Maxiter, sys_pmt, cmp_mesh, scft, scftAux, Mthd);
	
	dtol_energy = energy - old_energy;
	res = abs(dtol_energy);
	old_energy = energy;

	fprintf('\n &&&&&& ITER %d : Radius = %.15e, energy = %.15e, dtol_energy = %.15e, \n\n', ITER, Radius, energy, dtol_energy);
	fprintf('-----------------------------------------------------------------------------\n');
end
t2loop = clock;
time = etime(t2loop,t1loop); 

fprintf('========== END PROGRAM ==========\n\n\n');
fprintf('Time cost of SCFT ITERATION is " %.6f "  seconds\n\n\n', time);

close all;
figure(1)
showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));

saveName = sprintf('./results/phiA.%d', 000);
title(saveName);
set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
f = getframe(gcf);
imwrite(f.cdata,  [saveName, '.png']);
pause(1);

diary off;
