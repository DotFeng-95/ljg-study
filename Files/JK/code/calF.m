function F = calF(mymesh, mu, Sys_para)

chiN = Sys_para.Ndeg*Sys_para.chiAB;

mu1_int  = integrate_space(mymesh, mu(1,:));
mu2_int2 = integrate_w(mymesh, mu(2,:));

F = -mu1_int + mu2_int2/chiN;

end
