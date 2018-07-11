function Q = integrate_space(mymesh, q)
%%%%  how to cite the 'surface'

[lambda, wt] = quadpts(3);

nQuad = length(wt);
qq = q(mymesh.elem);

NT = size(mymesh.elem, 1);
b = zeros(NT, 1);
for i = 1:nQuad
    b = b + wt(i)*qq*lambda(i, :)';
end
b = b.*mymesh.area;
Q = sum(b);
