function [cmp_mesh, cmp_pmt, scft] = diblockmain(init_mesh, cmp_mesh, cmp_pmt, scft, scftAux, sys_pmt, Mthd)

nspecies = sys_pmt.Nspecies;
	chiN = sys_pmt.Ndeg*sys_pmt.chiAB;
	nt    = cmp_pmt.Nt;
	ndof  = cmp_pmt.Ndof;
	nelem = cmp_pmt.Nelem;
	tol   = cmp_pmt.TOL;
fprintf('========== Model Parameters ==========\n\n');
fprintf('chiN = %.3f, fA = %.3f\n\n', chiN, sys_pmt.fA);

fprintf('========== Discretization Points ==========\n\n');
fprintf('Dofs: [Nt, Nelem, Ndof] = [%d, %d, %d], PDEsolver: %s\n\n\n', nt, nelem, ndof, Mthd.PDE);

figure(1)
showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));
saveName = sprintf('./results/phiA.%d', 0);
title(saveName);
set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
f = getframe(gcf);
imwrite(f.cdata,  [saveName, '.png']);

dtol_energy = inf;
res = inf;

ITER = 0;

t1loop = clock;
[energy, scft] = evlSaddle(cmp_pmt, cmp_pmt.Maxiter, sys_pmt, cmp_mesh, scft, scftAux, Mthd);
old_energy = energy;

figure(2)
showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));
saveName = sprintf('./results/phiA.%d', 1);
title(saveName);
set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
f = getframe(gcf);
imwrite(f.cdata,  [saveName, '.png']);

fprintf('\n &&&&&& ITER %d : Radius = %.15e, energy = %.15e\n\n', ITER, cmp_mesh.radius, energy);

maxrefine = cmp_pmt.maxrefine;
refine = 1;
while true
    while res > cmp_pmt.TOL
        ITER = ITER + 1;

        radius = adaptSurface(cmp_mesh.radius, scft, cmp_pmt, init_mesh, Mthd, 1.0e-3);

        %%%%% generate computational meshgrid
        cmp_mesh = get_cmp_mesh(init_mesh, radius);
        [energy, scft] = evlSaddle(cmp_pmt, cmp_pmt.Maxiter, sys_pmt, cmp_mesh, scft, scftAux, Mthd);
        
        dtol_energy = energy - old_energy;
        res = abs(dtol_energy);
        old_energy = energy;

        fprintf('\n &&&&&& ITER %d : Radius = %.15e, energy = %.15e, dtol_energy = %.15e, \n\n', ITER, cmp_mesh.radius, energy, dtol_energy);
    end

    if refine > maxrefine
        break;
    else
        refine = refine + 1;
        [node, elem, ~, HB] = smeshuniformrefine(init_mesh.node, init_mesh.elem);
        node = init_mesh.surface.project(node);
        Ndof = cmp_pmt.Ndof;
        ndof = size(node, 1);
        mu = zeros(2, ndof);
        w = zeros(2, ndof);

        mu(1:2, 1:Ndof) = scft.mu;
        mu(1, Ndof+1:end) = (scft.mu(1, HB(:, 2)) + scft.mu(1, HB(:, 3)))/2;
        mu(2, Ndof+1:end) = (scft.mu(2, HB(:, 2)) + scft.mu(2, HB(:, 3)))/2;

        w(1,:) = mu(1,:) - mu(2,:);
        w(2,:) = mu(1,:) + mu(2,:);

        init_mesh.node = node;
        init_mesh.elem = elem;

        %%%%% generate computational meshgrid
        cmp_mesh = get_cmp_mesh(init_mesh, cmp_mesh.radius);
        cmp_pmt = get_cmp_pmt(sys_pmt, init_mesh);
        [scft, scftAux] = get_scft_data(sys_pmt, cmp_pmt);

        scft.w = w;
        scft.mu(1,:) = 0.5*(scft.w(1,:)+scft.w(2,:));
        scft.mu(2,:) = 0.5*(scft.w(2,:)-scft.w(1,:));

        scft.rho(1,:) = 0.5 + scft.mu(2,:)/chiN;
        scft.rho(2,:) = 1.0 - scft.rho(1,:);

        nspecies = sys_pmt.Nspecies;
        chiN = sys_pmt.Ndeg*sys_pmt.chiAB;
        nt    = cmp_pmt.Nt;
        ndof  = cmp_pmt.Ndof;
        nelem = cmp_pmt.Nelem;
        tol   = cmp_pmt.TOL;
        fprintf('========== Refine the init mesh ======\n\n');

        fprintf('========== Model Parameters ==========\n\n');
        fprintf('chiN = %.3f, fA = %.3f\n\n', chiN, sys_pmt.fA);

        fprintf('========== Discretization Points ==========\n\n');
        fprintf('Dofs: [Nt, Nelem, Ndof] = [%d, %d, %d], PDEsolver: %s\n\n\n', nt, nelem, ndof, Mthd.PDE);


        [energy, scft] = evlSaddle(cmp_pmt, cmp_pmt.Maxiter, sys_pmt, cmp_mesh, scft, scftAux, Mthd);
        old_energy = energy;

        dtol_energy = inf;
        res = inf;

        figure(2)
        showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));
        saveName = sprintf('./results/phiA.%d', 1);
        title(saveName);
        set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
        f = getframe(gcf);
        imwrite(f.cdata,  [saveName, '.png']);
    end
end



t2loop = clock;
time = etime(t2loop,t1loop); 

fprintf('========== END PROGRAM ==========\n\n\n');
fprintf('Time cost of SCFT ITERATION is " %.6f "  seconds\n\n\n', time);

figure(3)
showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));

saveName = sprintf('./results/phiA.%d', ITER);
title(saveName);
set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
f = getframe(gcf);
imwrite(f.cdata,  [saveName, '.png']);
pause(1);

end
