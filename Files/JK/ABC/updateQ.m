function sQ = updateQ(scft, cmpmesh,q_times_qplus)
%%%  compute single chain partition function sQ  

	nt = length(scft.qA(:,1));
	sQ(1) = integrate_space(cmpmesh, scft.q(nt,:));
	sQ(1) = sQ(1)/cmpmesh.tol_area;
    clear q_times_qplus;

end
