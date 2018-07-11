function sQ = updateQ(scft, cmpmesh, q_times_qplus)
%%%  compute single chain partition function sQ  

	nt = length(scft.q(:,1));
	sQ(1) = integrate_space(cmpmesh, scft.q(nt,:));
	sQ(1) = sQ(1)/cmpmesh.tol_area;

	%%%% TEST the code since Q is independent on s
%    Qarray  = zeros(1,nt);
%    for i =1:1:nt
%        Qarray(i) = integrate_space(cmpmesh, q_times_qplus(i,:));
%        Qarray(i) = Qarray(i)/cmpmesh.tol_area;
%    end
%    Qerr = max(Qarray)-min(Qarray);
%    figure(5);
%    plot(Qarray);

	clear q_times_qplus;

end
