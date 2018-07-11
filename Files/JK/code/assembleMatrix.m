function invMat = assembleMatrix(mymesh, Mthd)
%%%% all inverse matrice used in PDE solver of explicit schemes are calculated in advance
	[A, area] = stiff_matrix(mymesh.node, mymesh.elem);
	M = mass_matrix(mymesh.node, mymesh.elem, area);

	if strcmp(Mthd.PDE, 'ExplictEuler') 
		invMat = inv(M);
	elseif strcmp(Mthd.PDE, 'RK4')
		invMat = inv(M);
    elseif strcmp(Mthd.PDE, 'CN')
        invMat = [];
	elseif strcmp(Mthd.PDE, 'ImplicitEuler')
		invMat = [];
	end

	%    ML = 1./sum(M, 2);  /// gathering mass on diagonal line, second-order precsion
end
