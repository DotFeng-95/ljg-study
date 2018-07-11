function sys_pmt = get_sys_pmt()
%%%% set model parameters 

sys_pmt.Nspecies = 2;
sys_pmt.Nblend   = 1;
sys_pmt.Nblock   = 2;
sys_pmt.Ndeg  	 = 100;

%%%%  lamella phase
%sys_pmt.fA       = 0.5;
%sys_pmt.chiAB 	 = 0.15;

%%%%  cylinder phase
sys_pmt.fA       = 0.8;
sys_pmt.chiAB 	 = 0.30;

end
