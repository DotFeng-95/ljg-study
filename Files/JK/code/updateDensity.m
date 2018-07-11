function scft = updateDensity(scft, cmp_pmt, q_times_qplus)
%%%  compute densities

	index = 1;
	[index, scft.rho(1,:)] = Integration_t(index, cmp_pmt.Interval(1), q_times_qplus);
	[index, scft.rho(2,:)] = Integration_t(index, cmp_pmt.Interval(2), q_times_qplus);
	scft.rho = scft.rho/scft.sQ(1);

end
