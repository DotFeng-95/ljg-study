function scft = updatePropagator(cmp_pmt, cmpmesh, scft, Mthd)
	
	%%%  solve PDEs on surface 
	%%%  solve q_A 
	qA_index = 1;
	scft.qA(1,:) = 1.0;    
	[scft.qA,qA_index] = PDEsolver( cmp_pmt.Interval(1), cmpmesh, ... 
								   scft.w(1,:), scft.qA, ...
								   qA_index, Mthd);
    %%%  solve qAplus     
	qAplus_index = 1;
	scft.qAplus(1,:) = scft.qA(end,:);
	[scft.qAplus,qAplus_index] = PDEsolver( cmp_pmt.Interval(1), cmpmesh, ...
										   scft.w(1,:), scft.qAplus, ...
										   qAplus_index, Mthd );                          
                                         
    %%%  solve q_B                            
    qB_index = 1;
	scft.qB(1,:) = 1.0;
    [scft.qB,qB_index] = PDEsolver( cmp_pmt.Interval(2), cmpmesh, ...
								   scft.w(2,:), scft.qB, ...
								   qB_index, Mthd);
    %%%  solve qBplus
    qBplus_index = 1;
	scft.qBplus(1,:) = scft.qB(end,:);
	[scft.qBplus,qBplus_index] = PDEsolver( cmp_pmt.Interval(2), cmpmesh, .... 
										   scft.w(2,:), scft.qBplus, ...
										   qBplus_index, Mthd );                                
    %%%  solve q_C 
    qC_index = 1;
	scft.qC(1,:) = 1.0;
	[scft.qC,qC_index] = PDEsolver( cmp_pmt.Interval(3), cmpmesh, ...
								   scft.w(3,:), scft.qC, ...
								   qC_index, Mthd);
    %%%  solve qBplus
    qCplus_index = 1;
	scft.qCplus(1,:) = scft.qC(end,:);
	[scft.qCplus,qCplus_index] = PDEsolver( cmp_pmt.Interval(3), cmpmesh, .... 
										   scft.w(3,:), scft.qCplus, ...
										   qCplus_index, Mthd);

end
