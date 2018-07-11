function surfacedata = mcmullen

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'Hessian', @Hessian, ...
    'tanproperator', @tanproperator, 'initmesh',@initmesh, 'meancurvature', @meancurvature);



function [node,elem] = initmesh
M = load('meshdata/mcmullen.mat');
node = M.node;
elem = M.elem;
end

function z = phi(p)
% level set function
x = p(:,1);
y = p(:,2);
z = p(:,3);
x2 = x.^2;
y2 = y.^2;
z2 = z.^2;

z = (1 + x2).*(1 + y2).*(1 + z2) + 8*x.*y.*z -2;
end


function n = unitoutnormal(p)
% 
n = gradient(p);
l = sqrt(sum(n.^2,2));
n = n ./[l,l,l];

end

function n = gradient(p)
% 
x = p(:,1);
y = p(:,2);
z = p(:,3);
x2 = x.^2;
y2 = y.^2;
z2 = z.^2;
n = [2*x.*(1+y2).*(1+z2) + 8*y.*z,  2*y.*(1+x2).*(1+z2) + 8*x.*z, 2*z.*(1 + x2).*(1+y2) + 8*x.*y];

end

function [node,dist] = project(p)
% projection function 

s = sign(phi(p));

node = p;

normalAtNode = gradient(node);
valueAtNode = phi(node);
node = node - valueAtNode*ones(1,3).*normalAtNode./(dot(normalAtNode,normalAtNode,2)*ones(1,3));

vector = (-s)*ones(1,3).*(node - p);

d = s.*sqrt(dot(vector,vector,2));

normalAtNode = gradient(node);

node = p - d*ones(1,3).*normalAtNode./(sqrt(dot(normalAtNode,normalAtNode,2))*ones(1,3));

valueAtNode = phi(node);

normalAtNode = gradient(node);

vector = (-s)*ones(1,3).*(node - p);

d = s.*sqrt(dot(vector,vector,2));


e1 = normalAtNode./(sqrt(dot(normalAtNode, normalAtNode,2))*ones(1,3))-vector./(sqrt(dot(vector,vector,2))*ones(1,3));
error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));

k=1;
while max(abs(error)) > 1e-6 && k<200
    
    k=k+1;
    
    node = node - valueAtNode*ones(1,3).*normalAtNode./(dot(normalAtNode,normalAtNode,2)*ones(1,3));
    
    vector = -s*ones(1,3).*(node - p);
    d = s.*sqrt(dot(vector,vector,2));
    normalAtNode = gradient(node);
    
    node = p - d*ones(1,3).*normalAtNode./(sqrt(dot(normalAtNode,normalAtNode,2))*ones(1,3));
    
    valueAtNode = phi(node);
    
    normalAtNode = gradient(node);
    
    vector = (-s)*ones(1,3).*(node - p);
    
    d = s.*sqrt(dot(vector,vector,2));
    e1 = normalAtNode./(sqrt(dot(normalAtNode, normalAtNode,2))*ones(1,3))-vector./(sqrt(dot(vector,vector,2))*ones(1,3));
    error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));   
end


if nargout == 2
    dist = d;
end
end


function H = Hessian(p)

H = zeros(3,3,size(p,1));
end




function mc = meancurvature(p)
mc = 0;

end

function z = tanproperator(p)

z = 0;
end
end