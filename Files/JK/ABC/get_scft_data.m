function scft = get_scft_data(sys_pmt, cmp_pmt)

nspecies = sys_pmt.Nspecies;
nblend   = sys_pmt.Nblend;

ndof  = cmp_pmt.Ndof;

% number of dt of fA,fB,fC   
nt1 = cmp_pmt.Interval(1).n;
nt2 = cmp_pmt.Interval(2).n;
nt3 = cmp_pmt.Interval(3).n;

  % propagators equation
scft.qA     = zeros(nt1, ndof);
scft.qAplus = zeros(nt1, ndof);

scft.qB     = zeros(nt2, ndof);
scft.qBplus = zeros(nt2, ndof);

scft.qC     = zeros(nt3, ndof);
scft.qCplus = zeros(nt3, ndof);


scft.sQ    = zeros(1, nblend);
scft.rho   = zeros(nspecies, ndof);
scft.mu    = zeros(nspecies, ndof);

scft.w     = zeros(nspecies, ndof);
scft.grad  = zeros(nspecies, ndof);
end
