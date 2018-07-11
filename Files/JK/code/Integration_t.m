function [index, f] = Integration_t(index, interval, integrand)
%%%  fourth order quadrature scheme
%%%  \int f(t) dt = dt*(-(5/8)*(f_1+f_n) + (1/6)*(f_2+f_{n-1}) + \sum_{j=1}^n f_j)
%%%  INPUT: 
%%%      index: start integral point
%%%    integral: 
%%%   interval: integral interval, including dt and # dof of integral points

n  = interval.n;
dt = interval.dt;

f = -0.625* ( integrand(index,:) + integrand(index+n,:) ) + ...
	(1/6)* ( integrand(index+1,:) + integrand(index+n-1,:) ) - ...
	(1/24)* ( integrand(index+2,:) + integrand(index+n-2,:) );
f = f + sum(integrand(index:index+n,:), 1);
f = dt*f;

index = index+n;

end
