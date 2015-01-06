%% This function calcualtes the Jacobian of Bennett model. The differential is on a log scale
function [jac] = BennettJacobian(differential,K,t0,t,dt,L,x0)
% Get the original simulation
[T0,X0] =BennettModelSBML(K,t0,dt,0,x0);
x01 = X0(end,:);
[T,X] = BennettModelSBML(K,t,dt,L,x01);
traj_sim = X(:,6);
% initialize the matrices to calculate jacobian
leftJab = zeros(size(X,1),length(K)-3);
rightJab = zeros(size(X,1),length(K)-3);
% perform simulations
for i=1:size(leftJab,2)
    disp(strcat('Left i = ',num2str(i)));
    LeftParam = K;
    if i <=16
        LeftParam(i) = LeftParam(i)*10^(-differential);
        [T0,X0] =BennettModelSBML(LeftParam,t0,dt,0,x0);
        x01 = X0(end,:);
        [TL,XL] =BennettModelSBML(LeftParam,t,dt,L,x01);
        leftJab(:,i) = XL(:,6);
        
    else
        LeftParam(i+3) = LeftParam(i+3)*10^(-differential);
        [T0,X0] =BennettModelSBML(LeftParam,t0,dt,0,x0);
        x01 = X0(end,:);
        [TL,XL] =BennettModelSBML(LeftParam,t,dt,L,x01);
        leftJab(:,i) = XL(:,6);        
        
    end
    
end


for i =1:size(rightJab,2)
    disp(strcat('Right i = ',num2str(i)));
    RightParam = K;
    if i <=16
        RightParam(i) = RightParam(i)*10^(differential);
        [T0,X0] =BennettModelSBML(RightParam,t0,dt,0,x0);
        x01 = X0(end,:);
        [TR,XR] =BennettModelSBML(RightParam,t,dt,L,x01);
        rightJab(:,i) = XR(:,6);
        
    else
        RightParam(i+3) = RightParam(i+3)*10^(differential);
        [T0,X0] =BennettModelSBML(RightParam,t0,dt,0,x0);
        x01 = X0(end,:);
        [TR,XR] =BennettModelSBML(RightParam,t,dt,L,x01);
        rightJab(:,i) = XR(:,6);        
        
    end
end

jac = (rightJab-leftJab)/(2*differential);





end
