function area = spherearea(node, elem)

n1 = cross(node(elem(:, 3), :), node(elem(:, 2), :), 2);
n2 = cross(node(elem(:, 1), :), node(elem(:, 3), :), 2);
n3 = cross(node(elem(:, 2), :), node(elem(:, 1), :), 2);

l1 = sqrt(dot(n1, n1, 2));
l2 = sqrt(dot(n2, n2, 2));
l3 = sqrt(dot(n3, n3, 2));

n1 = n1./[l1, l1, l1];
n2 = n2./[l2, l2, l2];
n3 = n3./[l3, l3, l3];
a1 = pi - acos(dot(n2, n3, 2));
a2 = pi - acos(dot(n1, n3, 2));
a3 = pi - acos(dot(n1, n2, 2));

area = a1 + a2 + a3 - pi;
