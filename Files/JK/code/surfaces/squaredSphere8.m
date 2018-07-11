function surfacedata = squaredSphere8

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'Hessian', @Hessian, ...
    'tanproperator', @tanproperator, 'initmesh',@initmesh,...
    'meancurvature', @meancurvature, 'u',@u, 'gra_u',@gra_u, 'f',@f);



function [node,elem] = initmesh
M = load('meshdata/squaredsphere8_5847.mat');
node = M.node;
elem = M.elem;
end

function z = phi(p)
% level set function
x2 = p(:,1).^2;
y2 = p(:,2).^2;
z2 = p(:,3).^2;
x4 = x2.*x2;
y4 = y2.*y2;
z4 = z2.*z2;

z =  x4.*x4 + y4.*y4 + z4.*z4 - 1;
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
x4 = x2.*x2;
y4 = y2.*y2;
z4 = z2.*z2;
n = 8*[x4.*x2.*x, y4.*y2.*y, z4.*z2.*z];

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

    function z = u(p)
        p = project(p);
        x = p(:,1);y=p(:,2);z=p(:,3);
        z = sin(pi*x).*sin(pi*y).*sin(pi*z);
    end

    function z = f(p)
        p = project(p);
        x = p(:,1);y=p(:,2);z=p(:,3);
        x2 = x.^2; y2 = y.^2; z2 = z.^2;
        x4 = x2.^2; y4 = y2.^2; z4 = z2.^2;
        x6 = x4.*x2; y6 = y4.*y2; z6 = z4.*z2;
        x8 = x4.^2; y8 = y4.^2; z8 = z4.^2;
        x14 = x8.*x6; y14 = y8.*y6; z14 = z8.*z6;
        
        d = x14 + y14 + z14;
        d6 = z6 + y6; d8 = z8 + y8; d14 = z14 + y14;
        
        r = d6.*x14 + d14.*x6 + y6.*z6.*d8;
        r1 = 3.5 *r .* sin(pi*z) + pi*cos(pi*z).*z6.*z.*d;
        
        m1 = ((pi*d.^2.*sin(pi*z) + 3.5*r.*z6.*z.*cos(pi*z)).*sin(pi*y) + r1.*y6.*y.*cos(pi*y)).*sin(pi*x);
        m2 = (r1.*sin(pi*y) + pi*sin(pi*z).*y6.*y.*cos(pi*y).*d).*x6.*x.*cos(pi*x);
        z = 2*pi*(m1 + m2)./(d.^2);    
    end

    function z = gra_u(p)
        p = project(p);
        x = p(:,1);y=p(:,2);z=p(:,3);
        x2 = x.^2; y2 = y.^2; z2 = z.^2;
        x4 = x2.^2; y4 = y2.^2; z4 = z2.^2;
        x6 = x4.*x2; y6 = y4.*y2; z6 = z4.*z2;
        x7 = x6.*x; y7 = y6.*y; z7 = z6.*z;
        x14 = x7.^2; y14 = y7.^2; z14 = z7.^2;
        d0 = x14 + y14 + z14;
        
        r1 = pi*cos(pi*x).*sin(pi*y).*sin(pi.*z);
        r2 = pi*sin(pi*x).*cos(pi*y).*sin(pi.*z);
        r3 = pi*sin(pi*x).*sin(pi*y).*cos(pi.*z);
        
        d1 = (r1.*x7 + r2.*y7 + r3.*z7)./d0;
        
        z1 = r1 - d1.*x7;
        z2 = r2 - d1.*y7;
        z3 = r3 - d1.*z7;
        z = [z1, z2, z3];   
    end
        
end
