function [energy, scft] = evlSaddle(cmp_pmt, iterMax, sys_pmt, cmp_mesh, scft, scftAux, Mthd)

lid = fopen('./results/iteration.dat', 'a+');

tol = cmp_pmt.TOL;

res = inf;
Fold = inf;
ediff = inf;

iteration = 0;

while( (res > tol) && (iteration < iterMax) )

	t1loop = clock;
	iteration = iteration + 1;

	%%%  solve PDEs on surface 
	scft = updatePropagator(cmp_pmt, cmp_mesh, scft, Mthd);

    %%%  assemble integrand function for computing densities and Q
	q_times_qplus = scft.q.*flipdim(scft.qplus, 1);

    %%%  compute single chain partition function sQ  
	scft.sQ = updateQ(scft, cmp_mesh, q_times_qplus);
	
	%%%  compute densities
	scft = updateDensity(scft, cmp_pmt, q_times_qplus);

	%%%  update fields 
	[scft, err] = updateField(scft, Mthd.Iter, sys_pmt);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%   
 %%%  evaluate error
	err1 = err(1);  err2 = err(2);
    rho_err = max(scftAux.rho_old(:)-scft.rho(:));
    mu1_err  = max(scftAux.mu_old(1,:)-scft.mu(1,:));
    mu2_err  = max(scftAux.mu_old(2,:)-scft.mu(2,:));
    grad1_err  = max(scftAux.grad_old(1,:)-scft.grad(1,:));
    grad2_err  = max(scftAux.grad_old(2,:)-scft.grad(2,:));
    
	scftAux.rho_old = scft.rho;    
    scftAux.mu_old = scft.mu;
    scftAux.grad_old = scft.grad;

%    res = max(abs(err));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
	F = calF(cmp_mesh, scft.mu, sys_pmt);
	F = F/cmp_mesh.tol_area - log(scft.sQ(1));

	ediff = F-Fold;
	Fold = F;
	res = abs(ediff);
    
	fprintf('iteration %d:  sQ = %e, Energy = %.15e, ediff = %e, err1 = %e, err2 = %e\n', iteration, scft.sQ(1), F, ediff, err1, err2);
	fprintf(lid, '%d %.15e\t %.15e\t %.15e\t %.15e \t%.15e\n', iteration, scft.sQ(1), F, ediff, err1, err2);

	t2loop = clock;
	time = etime(t2loop,t1loop); 
	fprintf('Time cost is " %.6f "  seconds\n\n\n', time);

	if mod(iteration, 200)==0
		close all;
		figure(1)
		showresult(cmp_mesh.node, cmp_mesh.elem, scft.rho(1,:));

		saveName = sprintf('./results/phiA.%d', iteration);
		title(saveName);
		set(gcf, 'unit', 'normalized', 'position', [0.1, 0.1, 0.6, 0.25]);
		f = getframe(gcf);
		imwrite(f.cdata,  [saveName, '.png']);
		pause(1);
	end

end %% END while

energy = F;
fclose(lid);

end
