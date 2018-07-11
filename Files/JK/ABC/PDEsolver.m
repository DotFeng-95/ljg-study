function  [q, index] = PDEsolver(interval, mymesh, w, q, index, Mthd)
%%%%% explicit Euler method
	% INPUT: 
	%  SysParamts: 
	%   mymesh.node:
	%   mymesh.elem:
	%      w: Nx1
	%  index: point out the start position
	%
	%  OUTPUT:  	
	%      q: (nt+1)xN
	PDEflag = Mthd.PDE;
	nt = interval.n;
	dt = interval.dt;
	
	[A, area] = stiff_matrix(mymesh.node, mymesh.elem); 
	F = force_matrix(mymesh.node, mymesh.elem, w, area); 
	M = mass_matrix(mymesh.node, mymesh.elem, area);

	q0 = q(index,:);
    
	if strcmp(PDEflag, 'CN')
		AF = A + F;
        L = M + 0.5*dt*AF;
        for i = index:1:index+nt
            q(i, :) = q0;
            q1 = L\(-0.5*dt*AF*q0' + M*q0');
            q0 = q1';
        end
	end

	index = index+nt;
end
