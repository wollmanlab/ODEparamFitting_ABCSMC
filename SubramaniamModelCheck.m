%{
This function checks the validity of the Subramaniam Model by checking the conservations of the species.
It checks for the conservation principlles

Inputs:
    K = the vector of parameters
    x0 = initial conditions before ligand perturbation
    L = ligand amount for perturbation
    t0 = the time it takes for system to settle to steady state with zero input
    t = the time the system will go after t0
Outputs:
    T = vector of time points
    X = matrix of species concentrations

%}
function [] = SubramaniamModelCheck( K,x0,L,t0,t )
% The first plot panel has plots on conservation  
% Make subplot panel 
figure(1);
% Check for conservation of ligand
% Turn off the degradation flux of ligand
KL = K;
KL(44) = 0;
% first get the zero ligand input condition
[T0L,X0L,] = SubramaniamModel(KL,t0,x0);
xL = X0L(end,:);
xL(1) = L;
[TL,XL] = SubramaniamModel(KL,t,xL);
subplot(2,4,1);
plot(T0L,X0L(:,1) + X0L(:,3) + X0L(:,6) + X0L(:,8),'b',TL+ max(T0L),XL(:,1) + XL(:,3) + XL(:,6) + XL(:,8),'r');
title('Tot. Ligand');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% Check for conservation of receptor
% Turn off the degradation flux of receptor
KR = K;
KR(46) = 0;
[T0R,X0R] = SubramaniamModel(KR,t0,x0);
xR = X0R(end,:);
xR(1) = L;
[TR,XR] = SubramaniamModel(KR,t,xR);
subplot(2,4,2);
plot(T0R, (X0R(:,2)+ X0R(:,3) + X0R(:,6) + X0R(:,7) + X0R(:,8) + X0R(:,9) + X0R(:,10))./(X0R(1,2)+ X0R(1,3) + X0R(1,6) + X0R(1,7) + X0R(1,8) + X0R(1,9) + X0R(1,10)),'b',TR+ max(T0R),(XR(:,2)+ XR(:,3) + XR(:,6) + XR(:,7) + XR(:,8) + XR(:,9) + XR(:,10))./(X0R(1,2)+ X0R(1,3) + X0R(1,6) + X0R(1,7) + X0R(1,8) + X0R(1,9) + X0R(1,10)),'r');
title('Tot. Receptor');
xlabel('Time (seconds)');
ylabel('Fraction changed');
% first let the species go to steady states
% This steady state output can be used for all the subsequent 
[T0,X0,speciesArray] = SubramaniamModel(K,t0,x0);
% Check for conservation of G protein alpha
x1 = X0(end,:);
x1(1) = L;
[T1,X1,speciesArray1] = SubramaniamModel(K,t,x1);
subplot(2,4,3);
% use interpolation to get time series values in speciesArray that
% correspond to the values in the output X vector
GiDInterp0 = interp1(speciesArray(:,1),speciesArray(:,3),T0);
GiDInterp = interp1(speciesArray1(:,1), speciesArray1(:,3),T1);
% plot the conservation for G protien alpha
subplot(2,4,3);
plot(T0, (X0(:,11) + X0(:,12) + GiDInterp0)./(X0(1,11) + X0(1,12) + GiDInterp0(1)) ,'b',T1 + max(T0),(X1(:,11) + X1(:,12) + GiDInterp)./(X0(1,11) + X0(1,12) + GiDInterp0(1)),'r');
title('Tot. G alpha');
xlabel('Time (seconds)');
ylabel('Fraction changed');
% Check for conservation of G protein beta gamma
GRKGbyInterp0 = interp1(speciesArray(:,1),speciesArray(:,4),T0);
GRKGbyInterp = interp1(speciesArray1(:,1), speciesArray1(:,4),T1);
subplot(2,4,4);
plot(T0,(X0(:,4) + GiDInterp0 + GRKGbyInterp0)./(X0(1,4) + GiDInterp0(1) + GRKGbyInterp0(1)),'b', T1 + max(T0), (X1(:,4) + GiDInterp + GRKGbyInterp)./(X0(1,4) + GiDInterp0(1) + GRKGbyInterp0(1)),'r');
title('Tot. G beta gamma');
xlabel('Time (seconds)');
ylabel('Fraction changed)');
% Check for conservation of GRK
Cai2CaMGRKInterp0 = interp1(speciesArray(:,1),speciesArray(:,5),T0);
Cai2CaMGRKInterp  = interp1(speciesArray(:,1),speciesArray(:,5),T1);
subplot(2,4,5);
plot(T0,(X0(:,5) + GRKGbyInterp0 + Cai2CaMGRKInterp0)./(X0(1,5) + GRKGbyInterp0(1) + Cai2CaMGRKInterp0(1)),'b', T1 + max(T0), (X1(:,5) + GRKGbyInterp + Cai2CaMGRKInterp)./(X0(1,5) + GRKGbyInterp0(1) + Cai2CaMGRKInterp0(1)), 'r');
title('Tot. GRK');
xlabel('Time (seconds)');
ylabel('Fraction changed');
% Check for conservation of CaM
Cai2CaMInterp0 = interp1(speciesArray(:,1),speciesArray(:,6),T0); 
Cai2CaMInterp = interp1(speciesArray(:,1),speciesArray(:,6),T1);
subplot(2,4,6);
plot(T0, (X0(:,15) + Cai2CaMGRKInterp0 + Cai2CaMInterp0)./(X0(1,15) + Cai2CaMGRKInterp0(1) + Cai2CaMInterp0(1)),'b',T1 + max(T0), (X1(:,15) + Cai2CaMGRKInterp + Cai2CaMInterp)./(X0(1,15) + Cai2CaMGRKInterp0(1) + Cai2CaMInterp0(1)), 'r');
title('Tot. CaM');
xlabel('Time (seconds)');
ylabel('Fraction changed');
% Check for conservation of PIP2
ip3pInterp0 = interp1( speciesArray(:,1),speciesArray(:,2),T0); 
ip3pInterp = interp1( speciesArray(:,1),speciesArray(:,2), T1);
subplot(2,4,7);
plot(T0,(X0(:,13) + ip3pInterp0 + X0(:,14))./(X0(1,13) + ip3pInterp0(1) + X0(1,14)), 'b',T1 + max(T0),(X1(:,13) + ip3pInterp + X1(:,14))./(X0(1,13) + ip3pInterp0(1) + X0(1,14)), 'r');
title('Tot. PIP2');
xlabel('Time (seconds)');
ylabel('Fraction changed');
% The second plot panel has three plots 
figure(2);
% The second plot panel has plots on time profiles of species 
% Plot on ligand
subplot(3,3,1);
plot(T0,X0(:,1),'b',T1 + max(T0),X1(:,1),'r');
title('Ligand vs time');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% Plot on cytosolic calcium
subplot(3,3,2);
plot(T0,X0(:,16),'b',T1 + max(T0),X1(:,16),'r');
title('Cytosolic calcium vs time');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% Plot on ER calcium 
subplot(3,3,3);
plot(T0,X0(:,17),'b',T1 + max(T0),X1(:,17),'r');
title('ER calcium vs time');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% Plot on mitochondria calcium
subplot(3,3,4);
plot(T0,X0(:,19),'b',T1 + max(T0),X1(:,19),'r');
title('mitochondria calcium vs time');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% Plot on IP3 
subplot(3,3,5);
plot(T0,X0(:,19),'b',T1 + max(T0),X1(:,19),'r');
title('mitochondria calcium vs time');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% Plot on G protein beta gamma
subplot(3,3,6);
plot(T0,X0(:,4),'b',T1 + max(T0),X1(:,4),'r');
title('G protein beta gamma vs time');
xlabel('Time (seconds)');
ylabel('Concentration(micromolar)');
% The third plot panel has plots on fluxes

% Plot on IP3 generation flux

% Plot on Jch, the ip3 channel calcium flux

% Plot on Jserca

% Plot on Jer,leak

% Plot on Jpm,leak

% Plot on Jpm,ip3dep

% Plot on Jpmca

% Plot on Jncx

% Plot on Jmit,in 

% Plot on Jmit,out

end

