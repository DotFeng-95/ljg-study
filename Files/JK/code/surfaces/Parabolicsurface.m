function surfacedata = Parabolicsurface

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'Hessian', @Hessian, ...
    'tanproperator', @tanproperator, 'initmesh',@initmesh, 'meancurvature', @meancurvature);


function [node,elem] = initmesh(h)
[node,elem] = circlemesh(0,0,1,h);
x = node(:, 1);
y = node(:, 2);
z = x.^2 + y.^2;
newnode = [x, y, z];
node= project(newnode); 


function val = phi(p)
% level set function
val =  p(:,1).^ 2+ p(:,2).^ 2 - p(:, 3);


function n = unitoutnormal(p)
% 
s = 4.*p(:, 1).^2 + 4.*p(:, 2).^2+1
n = [2.*p(:, 1)./sqrt(s),2.*p(:, y)./sqrt(s), -1./sqrt(s)];

function n = gradient(p)
% 
x = p(:, 1);
y = p(:, 2);
z = p(:, 3);
n = [2*x, 2*y, -ones(length(x), 1)];

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

function H = Hessian(p)

H = zeros(3,3,size(p,1));
s = (4.*p(:, 1).^2 + 4.*p(:, 2).^2 + 1).^(1.5);
d11 = 2.*(4.*p(:, 2).^2 + 1)./s;
d12 = -8.*p(:, 1).*p(:, 2)./s;
d13 = 0;
d21 = -8.*p(:, 1).*p(:, 2)./s;
d22 = 2.*(4.*p(:, 1).^2 +1)./s;
d23 = 0
d31 = 4.*p(:, 1)./s;
d32 = 4.*p(:, 2)./s;
d33 = 0


for i=1:size(p,1)
    H(:,:,i)=[d11(i),d21(i),d31(i);d12(i),d22(i),d32(i);d13(i),d23(i),d33(i)];
end


function z = tanproperator(p)

z = zeros(3,3,size(p,1));
n = unitoutnormal(p);
for i=1:size(p,1)
   z(:,:,i)=n(i,:)'*n(i,:); 
end


