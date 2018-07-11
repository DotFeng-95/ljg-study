function w = initial_iield(w, mymesh, Sys_para, Cmp_para)
%%% initial fields 
%%% INPUT: 
%%%      mymesh: 

nspecies = Sys_para.Nspecies;
	ndof = Cmp_para.Ndof;

mu = zeros(nspecies, ndof);

%%%%  lamella phase
%%tmp = sin(10*pi*mymesh.node(:,1)).*sin(10*pi*mymesh.node(:,2)).*sin(10*pi*mymesh.node(:,3));
%tmp = sin(pi*mymesh.node(:,1));
%%tmp = ones(ndof, 1);
%mu(2,:) = tmp';

%%%%  cylinder phase
initsurf = spheresurface();
[inode,ielem] = initsurf.initmesh();
in = size(inode, 1);
mu(2,1:in) = 1.0;


w(1,:) = mu(1,:) - mu(2,:);
w(2,:) = mu(1,:) + mu(2,:);

clear mu;

end
