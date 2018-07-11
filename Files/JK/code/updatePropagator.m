function scft = updatePropagator(cmp_pmt, cmpmesh, scft, Mthd)
	
	%%%  solve PDEs on surface 
	%%%  solve q_A 
	q_index = 1;
	scft.q(1,:) = 1.0;

	[scft.q, q_index] = PDEsolver( cmp_pmt.Interval(1), cmpmesh, ... 
								   scft.w(1,:), scft.q, ...
								   q_index, Mthd);

    %%%  solve q_B 
	[scft.q, q_index] = PDEsolver( cmp_pmt.Interval(2), cmpmesh, ...
								   scft.w(2,:), scft.q, ...
								   q_index, Mthd);
    
	%%%  solve qplus_B
	qplus_index = 1;
	scft.qplus(1,:) = 1.0;
	[scft.qplus, qplus_index] = PDEsolver( cmp_pmt.Interval(2), cmpmesh, ...
										   scft.w(2,:), scft.qplus, ...,
										   qplus_index, Mthd );
	%%%  solve qplus_A
	[scft.qplus, qplus_index] = PDEsolver( cmp_pmt.Interval(1), cmpmesh, .... 
										   scft.w(1,:), scft.qplus, ...
										   qplus_index, Mthd );


end
