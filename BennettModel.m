function [ T,X ] = BennettModel(K,t,dt,L,x0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    % Rt, total no. of P2Y2 receptors
    Rt = K(1) ;
    % K1, unphosphorylated receptor dissociation constant
    K1 = K(2);
    % K2, phosphorylated receptor dissociation constant
    K2 = K(3);
    % kr, receptor recycling rate
    kr = K(4);
    %kp, receptor phosphorylation rate
    kp = K(5);
    % ke, receptor endocytosis rate
    ke = K(6);
    % xi, the fraction of mobile receptors 
    xi = K(7);
    % Gt, total number of G protein molecules
    Gt = K(8);
    % kdeg, IP3 degradation rate
    kdeg = K(9);
    % ka, G protein activation rate
    ka = K(10);
    % kd, G protein deactivation rate
    kd = K(11);
    % PIP2t, the total number of PIP2 molecules
    PIP2t = K(12);
    %rr, PIP2 replenishment rate
    rr = K(13);
    % delta, G protein intrinsic activity parameter
    delta = K(14);
    % Kc, the dissociation constant for calcium binding to PLC
    Kc = K(15);
    %alpha, effective signal gain parameter
    alpha = K(16);
    % v, cell volume 
    v = K(17);
    % Na, avogadro number
    Na = K(18);
    % epr, ratio of ER to cytosolic volume
    epr = K(19);
    % d1 , IP3 channel kinetic parameter 1
    d1 = K(20);
    % d2, IP3 channel kinetic parameter 2
    d2 = K(21);
    % d3, IP3 channel kinetic parameter 3
    d3 = K(22);
    % d5, IP3 channel kinetic parameter 5
    d5 = K(23);
    % a2, IP3 channel kinetic parameter 
    a2 = K(24);
    % Be, concentration of cytosolic endogeneous buffer
    Be = K(25);
    % Ke, cytosolic endogeneous buffer dissociation constant
    Ke = K(26);
    % Ber, concentration of ER endogeneous buffer
    Ber = K(27);
    % Ker, dissociation constant ER buffer
    Ker  = K(28);
    % Bx, concentration of cytosolic exogeneous buffer
    Bx = K(29);
    % Kx, cytosolic exogeneous buffer idssociation constant 
    Kx = K(30);
    % k3, pump dissociation constant
    k3 = K(31);
    % eta1, effective IP3 channel permeability
    eta1 = K(32);
    % eta2, effective ER leak permeability
    eta2 = K(33) ;
    % eta3, effective calcium pump permeability
    eta3 = K(34);
    % Cat, total concentration of calcium
    Cat = K(35);
    Tgrid = (0:dt:t)';
    function [dx_dt] = F(t,x)
        dx_dt = zeros(7,size(x,2));
        %{
        % x0(1) = initial concentration for unphosphorylated  surface receptors
x0(1) = K(1)*K(7);
% x0(2) = initial concentration for phosphorylated surface receptors
x0(2) = 0;
% x0(3) = initial concentration of activated G protein
x0(3) = K(36);
%x0(4) = PIP2 initial concentration
x0(4) = K(37);
% x0(5) = IP3 initial concentration
x0(5) = K(38);
% x0(6) = cytosolic calcuim initial concentration
x0(6) = K(39);
% initial concentration of h
x0(7) = K(40); 
  %}    
        
        % rho_r is the ratio of teh number of ligand bound receptors to the
        % total number of receptors
        rho_r = L*x(1,:)./(xi*Rt*(K1 + L));
        % rh is the rate coefficient for the hydrolysis of PIP2
        rh = alpha*(x(6,:)./(Kc + x(6,:))).*x(3,:);
        % Beta is the buffering function
        Beta = 1/(1 + Ke*Be/(Ke + x(6,:)).^2 + Kx*Bx./(Kx + x(6,:)).^2);
        % zeta
        zeta = d2.*(x(5,:) + d1)./(x(5,:) + d3);
        % tauh is the time constant for IP3 receptor inactivation 
        tauh = 1./(a2.*(zeta + x(6,:))); 
        %hinf is used in the equaiton for the rate of change of h
        hinf = zeta./(zeta + x(6,:));
        %minf is used in the cytosolic calcium equation
        minf = (x(5,:)./(d1 + x(5,:))).*(x(6,:)./(d5 + x(6,:)));
        % gamma is used in the equation to calculate Caer
        gamma = 1./(1 + Be./(Ke + x(6,:)) + Bx./(Kx + x(6,:)));
        % Caer is the er calcium
        Caer = Ker.*(Cat - x(6,:)./gamma)./(Ber*epr);        
        % rate of change of unphosphorylated surface receptors
        dx_dt(1,:) = kr*Rt - (kr + kp*L/(K1 + L)).*x(1,:) - kr*x(2,:);
        % rate of change of phosphorylated surface receptors
        dx_dt(2,:) = L*(kp*x(1,:)./(K1 + L) - ke*x(2,:)./(K2 + L));
        % rate of change of activated G protein 
        dx_dt(3,:) = ka*(delta + rho_r).*(Gt - x(3,:)) - kd*x(3,:);
        % rate of change of PIP2 
        dx_dt(4,:) = -(rh + rr).*x(4,:) - rr*Na*v*x(5,:) + rr*PIP2t;
        % rate of change of IP3 
        dx_dt(5,:) = rh*x(4,:)./(Na*v) -kdeg*x(5,:);
        % rate of change of cytosolic calcium
        dx_dt(6,:) = Beta*(epr*(eta1.*minf.^3.*x(7,:).^3 + eta2).*(Caer - x(6,:)) - eta3.*(x(6,:).^2./(k3^2 + x(6,:).^2)));
        % rate of change of h
        dx_dt(7,:) = (hinf -x(7,:))./tauh;
        
    
    end
    options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Refine',1000);
    [T,X] = ode15s(@F,Tgrid,x0,options);

end

