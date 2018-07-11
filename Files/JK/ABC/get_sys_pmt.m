function sys_pmt = get_sys_pmt()
%%%% set ABC star model parameters 

sys_pmt.Nspecies = 3;
sys_pmt.Nblend   = 1;
sys_pmt.Nblock   = 2;
sys_pmt.Ndeg  	 = 100;


%%%%  parameters
sys_pmt.fA       = 0.2;
sys_pmt.fB       = 0.4;
sys_pmt.fC       = 0.4;

sys_pmt.chiAB    = 0.11;
sys_pmt.chiBC    = 0.25;
sys_pmt.chiAB 	 = 0.30;

end
