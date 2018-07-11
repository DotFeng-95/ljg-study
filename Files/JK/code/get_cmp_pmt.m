function cmp_pmt = get_cmp_pmt(sys_pmt, cmp_mesh)
%%%% get the computational parameters

          cmp_pmt.dim   = 2;
		 cmp_pmt.Ndof   = size(cmp_mesh.node, 1);
		 cmp_pmt.Nelem  = size(cmp_mesh.elem, 1);
	    cmp_pmt.dtMax   = 5.0e-3;
cmp_pmt.Interval(1).a   = 0.0;
cmp_pmt.Interval(1).b   = sys_pmt.fA;
%%%%  cmp_pmt.Interval(1).n: record the number of segments 
cmp_pmt.Interval(1).n   = ceil(sys_pmt.fA / cmp_pmt.dtMax);
cmp_pmt.Interval(1).dt  = (cmp_pmt.Interval(1).b - cmp_pmt.Interval(1).a) / cmp_pmt.Interval(1).n;
cmp_pmt.Interval(2).a   = sys_pmt.fA;
cmp_pmt.Interval(2).b   = 1.0;
%%%%  cmp_pmt.Interval(2).n: record the number of segments 
cmp_pmt.Interval(2).n   = ceil((1.0-sys_pmt.fA) / cmp_pmt.dtMax);
cmp_pmt.Interval(2).dt  = (cmp_pmt.Interval(2).b - cmp_pmt.Interval(2).a) / cmp_pmt.Interval(2).n;
		    cmp_pmt.Nt  = cmp_pmt.Interval(1).n +  cmp_pmt.Interval(2).n + 1;
		    cmp_pmt.TOL = 1.0e-6; %% 
		cmp_pmt.Maxiter = 5000;
	 	 cmp_pmt.Initer = 200;

end
