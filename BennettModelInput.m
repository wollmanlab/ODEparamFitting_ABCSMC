%% This function generates the input parameters and ICs for the Bennett model
% Units of concentration are in micromolar, whereas units of time are in
% seconds, lengths are in micron
function [K,t0,t,dt,x0] = BennetModelInput()
    % initialize the parameter for 
    K = zeros(35,1);

    % Rt, total no. of P2Y2 receptors
    K(1) = 2e4;
    % K1, unphosphorylated receptor dissociation constant
    K(2) = 5;
    % K2, phosphorylated receptor dissociation constant
    K(3) = 100;
    % kr, receptor recycling rate
    K(4) = 1.75e-4;
    %kp, receptor phosphorylation rate
    K(5) = 0.03;
    % ke, receptor endocytosis rate
    K(6) = 6e-3;
    % xi, the fraction of mobile receptors 
    K(7) = 0.85;
    % Gt, total number of G protein molecules
    K(8) = 1e5;
    % kdeg, IP3 degradation rate
    K(9) = 1.25;
    % ka, G protein activation rate
    K(10) = 0.017;
    % kd, G protein deactivation rate
    K(11) = 0.15;
    % PIP2t, the total number of PIP2 molecules
    K(12) = 5e4;
    %rr, PIP2 replenishment rate
    K(13) = 10;
    % delta, G protein intrinsic activity parameter
    K(14) = 1.234e-3;
    % Kc, the dissociation constant for calcium binding to PLC
    K(15) = 0.4;
    %alpha, effective signal gain parameter
    K(16) = 2.781e-5;
    % v, cell volume. The paper lists this constant in unit of meter. Here it is listed in unit of liter 
    K(17) = 5e-13;
    % Na, avogadro num
    K(18) = 6.02252e17;
    % epr, ratio of ER to cytosolic volume
    K(19) = 0.185;
    % d1 , IP3 channel kinetic parameter 1
    K(20) = 0.13;
    % d2, IP3 channel kinetic parameter 2
    K(21) = 1.05;
    % d3, IP3 channel kinetic parameter 3
    K(22) = 0.943;
    % d5, IP3 channel kinetic parameter 5
    K(23) = 0.0823;
    % a2, IP3 channel kinetic parameter 
    K(24) = 0.2;
    % Be, concentration of cytosolic endogeneous buffer
    K(25) = 150;
    % Ke, cytosolic endogeneous buffer dissociation constant
    K(26) =  10;
    % Ber, concentration of ER endogeneous buffer
    K(27) = 120000;
    % Ker, dissociation constant ER buffer
    K(28) = 1200;
    % Bx, concentration of cytosolic exogeneous buffer
    K(29) = 50;
    % Kx, cytosolic exogeneous buffer idssociation constant 
    K(30) = 0.2;
    % k3, pump dissociation constant
    K(31) = 0.4;
    % eta1, effective IP3 channel permeability
    K(32) = 575;
    % eta2, effective ER leak permeability
    K(33) = 5.2;
    % eta3, effective calcium pump permeability
    K(34) = 45;
    % Cat, total concentration of calcium
    K(35) = 67;
    
    % Cabas, basal calcium concentration
    Cabas = K(31)/(K(34)*K(27)/(K(33)*K(28)*K(35)) -1)^(1/2);
    % Gbas, basal number of G protein molecules
    Gbas = K(10)*K(14)*K(8)/(K(10)*K(14) + K(11));
    %rhbas is the basal level of rh
    rhbas = K(16)*(Cabas/(K(15) + Cabas))*Gbas;
    % PIP2bas, basal number of PIP2 molecules
    PIP2bas = K(5)*K(13)*K(12)/(K(5)*(rhbas + K(13)) + rhbas*K(13));
    % IP3bas, basal IP3 concentration
    IP3bas = rhbas*K(13)*K(12)/(K(17)*K(18)*(K(5)*rhbas + K(13)) + rhbas*K(13));
    % hbas
    hbas =K(21)*(K(20) + IP3bas)/(K(20)*K(21) + IP3bas*(Cabas + K(21)) + K(22)*Cabas);

    x0 = zeros(7,1);
    % x0(1) = initial concentration for unphosphorylated  surface receptors
    x0(1) = K(7)*K(1);
    % x0(2) = initial concentration for phosphorylated surface receptors
    x0(2) = 0;
    % x0(3) = initial concentration of activated G protein
    x0(3) = Gbas;
    %x0(4) = PIP2 initial concentration
    x0(4) =PIP2bas;
    % x0(5) = IP3 initial concentration
    x0(5) =IP3bas;
    % x0(6) = cytosolic calcuim initial concentration
    x0(6) =Cabas;
    % initial concentration of h
    x0(7) = hbas;






    % t0 is the time the system runs before perturbation
    t0 = 0;
    %
    t = 240;
    %
    dt = 1;


end