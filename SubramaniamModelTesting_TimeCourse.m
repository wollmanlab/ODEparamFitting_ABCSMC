%% This script reproduces the time course data in Figure 2 in the Subramaniam paper (2007) 

%% Create model inputs (using the Subramniam Input function) 
% Inputs:
%     K = the vector of parameters
%     x0 = initial conditions before ligand perturbation
%     L = ligand amount for perturbation
%     t0 = the time it takes for system to settle to steady state with zero input
%     t = the time the system will go after t0
% Outputs:
%     T = vector of time points
%     X = matrix of species concentrations

%% Figure #1 Receptor time course data
% create inputs
[K,x0,L,t0,t,dt] = SubramaniamInput; 
% first get the zero ligand input condition
[T0L,X0L,speciesArray0,fluxMatrix0] = SubramaniamModel(K,t0,x0,dt);
xL = X0L(end,:);
xL(1) = L;
[TL,XL,speciesArray1,fluxaMatrix1] = SubramaniamModel(K,t,xL,dt);

% plot the R/Rtot fraction
figure(1);
clf;
subplot(2,2,1);
% Rtot0 is the total surface receptors before perturbation
Rtot0 = X0L(:,2) + X0L(:,3) + X0L(:,6) + X0L(:,7) + X0L(:,8) + X0L(:,9) + X0L(:,10);
%RtotL is the total surface receptors after perturbation
RtotL = XL(:,2) + XL(:,3) + XL(:,6) + XL(:,7) + XL(:,8) + XL(:,9) + XL(:,10);
plot(T0L,X0L(:,2)./Rtot0,'b', TL + max(T0L),XL(:,2)./RtotL,'r');
title('R/Rtot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');

%plot the L.R/Rtot fraction
subplot(2,2,2);
plot(T0L,X0L(:,3)./Rtot0,'b', TL + max(T0L),XL(:,3)./RtotL,'r');
title('[L.R]/Rtot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');

% plot the [L.Ri]/Rtot fraction
subplot(2,2,3);
plot(T0L,X0L(:,8)./Rtot0,'b', TL + max(T0L),XL(:,8)./RtotL,'r');
title('[L.Ri]/Rtot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');

% plot the [Rpi]/Rtot fraction 
subplot(2,2,4);
plot(T0L,X0L(:,9)./Rtot0,'b', TL + max(T0L),XL(:,9)./RtotL,'r');
title('[Rpi]/Rtot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');

%% Fig #2 G protein complexes
figure(2);
clf;
subplot(2,3,1);
% total G protein related complexes before perturbation
Gbytot0 = X0L(:,4) + speciesArray0(:,2) + speciesArray0(:,3);
Gaitot0 = X0L(:,12) + X0L(:,11) + speciesArray0(:,2);
Gbytot1 = XL(:,4) + speciesArray1(:,2) + speciesArray1(:,3);
Gaitot1 = XL(:,12) + XL(:,11) + speciesArray1(:,2);
% plot the ratio of Gby with Gbytot
plot(T0L,X0L(:,4)./Gbytot0,'b',TL + max(T0L),XL(:,4)./Gbytot1,'r');
title('[Gby]/Gbytot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');
% plot the ratio of GaiT with Gaitot
subplot(2,3,2);
plot(T0L,X0L(:,11)./Gaitot0,'b',TL + max(T0L),XL(:,11)./Gaitot1,'r');
title('[GaiT]/Gaitot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');
% plot the ratio of GRK.Gby with Gbytot
subplot(2,3,3);
plot(T0L,X0L(:,3)./Gbytot0,'b',TL + max(T0L),XL(:,3)./Gbytot1,'r');
title('[GRK.Gby]/Gbytot');
axis([0 max(TL)+max(T0L) 0 1]);
xlabel('Time (seconds)');
ylabel('Fraction');

%% Fig #3 Q and h
figure(3);
clf;
% plot Q
subplot(1,2,1);
plot(T0L,speciesArray0(:,6),'b',TL + max(T0L),speciesArray1(:,6),'r');
title('Q');
xlabel('Time (seconds)');
% plot h
subplot(1,2,2);
plot(T0L,X0L(:,18),'b',TL + max(T0L),XL(:,18),'r');
title('h');
xlabel('Time (seconds)');
%% Fig #4 IP3 and PIP2
figure(4);
clf;
% plot ip3
subplot(1,2,1);
plot(T0L,X0L(:,14),'b',TL+ max(T0L),XL(:,14),'r');
title('IP3');
xlabel('Time (seconds)');
% plot PIP2/PIP2tot fraction
subplot(1,2,2);
%PIP2tot0 is the PIP2tot before perturbation  
PIP2tot0 = speciesArray0(:,1) + X0L(:,13) + X0L(:,14);
%PIP2tot1 is the PIP2tot after perturbation
PIP2tot1 = speciesArray1(:,1) + XL(:,13) + XL(:,14);
plot(T0L,X0L(:,13)./PIP2tot0,'b',TL + max(T0L),XL(:,13)./PIP2tot1);
title('PIP2/PIP2tot');
xlabel('Time (seconds)');
ylabel('Fraction');
axis([0 max(TL)+max(T0L) 0 1]);

%% Fig #5 cytosolic Ca, Caer and Camit
figure(5);
clf;
%plot Cai
subplot(1,3,1);
plot(T0L,X0L(:,16),'b',TL+ max(T0L),XL(:,16),'r');
title('cytosolic Ca');
xlabel('Time (seconds)');
ylabel('Concentration (micromolar)');

% plot Caer
subplot(1,3,2);
plot(T0L,X0L(:,17),'b',TL + max(T0L),XL(:,17),'r');
title('ER Ca');
xlabel('Time (seconds)');
ylabel('Concentration (micromolar)');

%plot Ca mit
subplot(1,3,3);
plot(T0L,X0L(:,19),'b',TL + max(T0L),XL(:,19),'r');
title('Mit Ca');
xlabel('Time (seconds)');
ylabel('Concentration (micromolar)');
%% Fig #6 plot the evolution of G protein complexes and their conservation 
figure(6);
clf;


%