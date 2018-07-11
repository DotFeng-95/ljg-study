function n = countdots(node, elem, rho)

[m1, m2] = size(rho);
if m1 < m2
    rho = rho';
end

rho1 = (rho(elem(:, 1)) + rho(elem(:, 2)) + rho(elem(:, 3)))/3;
elem1 = elem(rho1 > 0.5, :);

T = auxstructure(elem1);
edge = T.edge;
clear T;

F = size(elem1, 1);
E = size(edge, 1);
V = length(unique(elem1(:)));
n = V + F - E;
end

