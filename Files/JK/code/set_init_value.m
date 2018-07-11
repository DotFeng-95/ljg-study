function u = set_init_value(node, elem, n)

N = size(node, 1);
T = auxstructure(elem);
edge = double(T.edge);
clear T

p2p = sparse(edge, edge(:, [2,1]), 1, N, N);

u = zeros(N, 1);
u(1:12) = 1;

for i=1:n
    isOne = p2p*u > 0;
    u(isOne & u == 0) = 1;
end
end
