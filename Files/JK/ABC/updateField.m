function [scft, err] = updateField(scft, flag, sys_para)


%% give rho initial value to compute mu1 and mu2

scft.rho(1,:) = 0.2
scft.rho(2,:) = 0.3
scft.rho(3,:) = 0.4

chiAB = sys_para.Ndeg*sys_para.chiAB;
chiAC = sys_para.Ndeg*sys_para.chiAC;
chiBC = sys_para.Ndeg*sys_para.chiBC;
Ndeg = sys_para.Ndeg;
alpha = (chiAC+chiAB-chiBC)/2*chiAC

delta = chiAB^2 + chiAC^2 + chiBC^2 - 2*chiAB*chiAC -2*chiAB* chiBC-2*chiAC*chiBC;
Xi1 = -delta/(4*chiAC);
Xi2 = chiAC;

sigm1A = 1/3;sigm1B = -2/3;sigm1C=1/3;
sigm2A = (1+ alpha)/3; sigm2B = (1-2*alpha)/3; sigma2C= (alpha-2)/3;


scft.mu(2,:) = 2*Ndeg*Xi1*(sigm1A*scft.rho(1,:)+sigm1B*scft.rho(2,:)+sigm1C*scft.rho(3,:));
scft.mu(3,:) = 2*Ndeg*Xi2*(sigm2A*scft.rho(1,:)+sigm2B*scft.rho(2,:)+sigm2C*scft.rho(3,:));

scft.w(1,:) = scft.mu(1,:) -sigm1A*scft.mu(2,:)%% 未完成







	if strcmp(flag, 'Euler')

%%% step length for update fields
		lambda = [2, 2, 2];
%%%  Euler scheme for updating fields mu

		scft.grad(1,:) = scft.rho(1,:)+scft.rho(2,:)+ scft.rho(3,:)-1.0;
        scft.grad(2,:) = scft.mu(2,:)/(2*Ndeg*simge1) - 1/3*scft.rho(1,:) + 2/3 * scft.rho(2,:)-1/3*scft.rho(3,:); 
        scft.grad(3,:) = scft.mu(2,:)/(2*Ndeg*simge1) - (1+alpha)/3*scft.rho(1,:) + (1-2*alpha)/3*scft.rho(2,:) - (alpha-2)/3*scft.rho(3,:); 
		err(1) = max(norm(scft.grad(1,:), inf));
		err(2) = max(norm(scft.grad(2,:), inf));
        err(3) = max(norm(scft.grad(3,:), inf));
        
    
        scft.w(1,:) = scft.w(1,:) + lambda(1)*scft.grad(1,:);
        scft.w(2,:) = scft.w(2,:) + lambda(2)*scft.grad(2,:);
        scft.w(3,:) = scft.w(3,:) + lambda(3)*scft.grad(3,:);
        
        
		scft.w(1,:) = scft.mu(1,:) - scft.mu(2,:);
		scft.w(2,:) = scft.mu(1,:) + scft.mu(2,:);
	end
end