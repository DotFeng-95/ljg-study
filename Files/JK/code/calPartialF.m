function partialF = calPartialF(radius, scft, cmp_pmt, Init_mesh, Mthd)	

	partialF = 0;
	cmpmesh = get_cmp_mesh(Init_mesh, radius);

	scft = updatePropagator(cmp_pmt, cmpmesh, scft, Mthd);
	q_times_qplus = scft.q.*flipdim(scft.qplus, 1);

	scft.sQ = updateQ(scft, cmpmesh, q_times_qplus);

	partialF = partialF - sum(log(scft.sQ(:)));

end
