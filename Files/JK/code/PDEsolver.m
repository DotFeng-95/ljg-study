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

	N = size(mymesh.node, 1);
	nt = interval.n;
	dt = interval.dt;
	
	[A, area] = stiff_matrix(mymesh.node, mymesh.elem); 
	F = force_matrix(mymesh.node, mymesh.elem, w, area);
	M = mass_matrix(mymesh.node, mymesh.elem, area);

	option.solver = 'CG';
	option.tol = 1.0e-6;
	option.printlevel = 0;

	q0 = q(index,:);

	if strcmp(PDEflag, 'ExplictEuler')
		AF = A+F;
		for i=index:1:index+nt   
			q(i,:) = q0;
			b = -dt*AF*q0';
%            q1 = Mthd.invMat*b + q0';
			q1 = M\b + q0';
			q0 = q1';
%            showresult(mymesh.node, mymesh.elem, q0)
		end
		
	elseif strcmp(PDEflag, 'ImplicitEuler')
		MAF = M+dt*(A+F); 
	%    MEL = inv(M+dt*(A+F));

		for i=index:1:index+nt   
			q(i,:) = q0;
%%%%%   提前算好右端矩阵的逆
%            q1 = MEL*M*q0';
%%%%%   使用代数多重网格求线性代数方程组
%            [q1, info] = amg(MAF, M*q0', option);                 
			q1 = MAF\(M*q0');                 
			q0 = q1';
%            showresult(mymesh.node, mymesh.elem, q0)
%            pause();
		end
    elseif strcmp(PDEflag, 'CN')
        AF = A + F;
        CM = sum(M, 2);
        for i=index:1:index+nt
            q(i, :) = q0;
            q1 = q0;
            err = inf;
            while err > 1e-8
                b = -dt*AF*(q1'+ q0')/2;
                q2 = b./CM +q0';
                err = norm(q2'-q1);
                q1 = q2';
            end
            q0 = q1;
        end
    elseif strcmp(PDEflag, 'ImplicitCN')
        AF = A + F;
        L = M + 0.5*dt*AF;
        for i = index:1:index+nt
            q(i, :) = q0;
            q1 = L\(-0.5*dt*AF*q0' + M*q0');
            q0 = q1';
        end
    elseif strcmp(PDEflag, 'ImplicitCN1')
        AF = A + F;
        CM = sum(M, 2);
        N = length(CM);
        L = spdiags(CM, 0, N, N) + 0.5*dt*AF;
        for i = index:1:index+nt
            q(i, :) = q0;
            q1 = L\(-0.5*dt*AF*q0' + M*q0');
            q0 = q1';
        end
	elseif strcmp(PDEflag, 'RK4')
		B = Mthd.invMat*(A+F);
		for i=index:1:index+nt
			q(i,:) = q0;

			k1  = q0;
			k2  = q0 + 0.5*dt*k1;
			k3  = q0 + 0.5*dt*k2;
			k4  = q0 + dt*k3;
			q1 = q0' + (dt/6)*B* (k1+ 2*k2 + 2*k3 + k4)';

			q0  = q1';
		end
	end

	index = index+nt;
end
