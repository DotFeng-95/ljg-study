function W = integrate_w(mymesh, w)

[lambda, wt] = quadpts(3);

nQuad = length(wt);
ww = w(mymesh.elem);

NT = size(mymesh.elem, 1);
b = zeros(NT, 1);
for i = 1:nQuad
    b = b + wt(i)*(ww*lambda(i, :)').^2;
end
b = b.*mymesh.area;
W = sum(b);
