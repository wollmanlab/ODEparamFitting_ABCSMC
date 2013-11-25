%% This script tests the accuracy of the Bennett Model by reproducing the plots in the paper
% 
ver =1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
L = 0;
%x0 = [2e4*0.85; 0; 0 ; 5e4 ; 0; 0.05; 1];
t0 = 3000;
t = 2100;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
% plot the variables to make sure that they run to equilibrium
figure(1);
subplot(2,4,1);
plot(T0,X0(:,1));
title('Unphosphorylated receptors');
xlabel('Time (seconds)');
ylabel('number');

subplot(2,4,2);
plot(T0,X0(:,2));
title('Phorylated receptors');
xlabel('Time (seconds)');
ylabel('number');

subplot(2,4,3);
plot(T0,X0(:,3));
title('activated G protein');
xlabel('Time (seconds)');
ylabel('number');

subplot(2,4,4);
plot(T0,X0(:,4));
title('PIP2');
xlabel('Time (seconds)');
ylabel('number');

subplot(2,4,5);
plot(T0,X0(:,5));
title('IP3');
xlabel('Time (seconds)');
ylabel('micromolar');

subplot(2,4,6);
plot(T0,X0(:,6));
title('Cytosolic calcium');
xlabel('Time (seconds)');
ylabel('micromlar');

subplot(2,4,1);
plot(T0,X0(:,7));
title('h');
xlabel('Time (seconds)');
ylabel('fraciton');


%% Reproduce Fig.3 (a) in the paper that relates ligand concentration to the percentage of surface receptors at equilibrium
ver = 1;
receptorArray = zeros(5,1); 
utpArray = [10^-1; 1; 10; 10^2 ; 10^3];
[K,t0,t,dt,x0] = BennettModelInput(ver);
t0 = 3000;
t = 2100;
L = 0;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
%x0(5) = 0.01;
% Get data points for the receptor array to be plotted
for i = 1:size(receptorArray,1)
    [T,X] = BennettModel(K,t,dt,utpArray(i),x01);
    receptorArray(i) = (X(end,1) + X(end,2))/K(1);
end
figure(2);
clf;
subplot(1,2,1);
plot(log10(utpArray),receptorArray);
title('Surface Receptors at Equilibrium');
xlabel('log 10 ligand (micromolar)');
ylabel('Surface Receptors Percentage');


L = 10^3;
[T,X] = BennettModel(K,t,dt,L,x01);
subplot(1,2,2);
surfRecep = (X(:,1) + X(:,2))/K(1);
plot(T,surfRecep);
title('Surface receptor time course');
xlabel('Time (second)');
ylabel('Surface receptor( fraction)');


%% Reproduce Fig. 5 in the Bennett paper
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
%x0(5) = 0.01;
L = 0;
t0 = 3000;
t = 480;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
rrArray = [10;0.1;0.015];
tArray = cell(3,1);
xArray = cell(3,1);
% Perform simulations
for i = 1:size(rrArray,1)
    K1 = K;
    K1(13) = rrArray(i);
    L = 10^3;
    [T,X] = BennettModel(K1,t,dt,L,x01);
    tArray{i} = T; 
    xArray{i} = X;
end
% Plot the figures
figure(3);
% plot Figure 5 (a)
subplot(2,2,1);
xArray1 = xArray{1};
xArray2 = xArray{2};
xArray3 = xArray{3};
tArray1 = tArray{1};
tArray2 = tArray{2};
tArray3 = tArray{3};
plot(tArray1(1:end/2),xArray1(1:end/2,5),'r',tArray2(1:end/2),xArray2(1:end/2,5),'b',tArray3(1:end/2),xArray3(1:end/2,5),'g');
title( 'IP3 variaiton with rr');
xlabel('Time (seconds)');
ylabel('IP3 (micromolar)');
legend('rr = 10','rr = 0.1','rr = 0.015');
% plot Figure 5(b)
subplot(2,2,2);
plot(tArray1(1:end/2),xArray1(1:end/2,6),'r',tArray2(1:end/2),xArray2(1:end/2,6),'b',tArray3(1:end/2),xArray3(1:end/2,6),'g');
title(' Calcium variation with rr');
xlabel('Time (seconds)');
ylabel('Calcium (micromolar)');
legend('rr = 10','rr = 0.1', 'rr = 0.015');
% plot figure 5(c)
subplot(2,2,3);
plot(tArray1,xArray1(:,4)/K(12),'r',tArray2,xArray2(:,4)/K(12),'b',tArray3,xArray3(:,4)/K(12),'g');
title(' PIP2 variation with rr');
xlabel('Time (seconds)');
ylabel('Relative PIP2 fraction');
legend('rr = 10','rr = 0.1', 'rr = 0.015');
% plot figure 5 (d)
subplot(2,2,4);
plot(tArray1(1:end/2),xArray1(1:end/2,3)/K(8));
title('Fraction of activated G protein');
xlabel('Time (seconds)');
ylabel('Relative G protein fraction');
%% This code block plots the time series of receptor variables and their conservation
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
%x0 = [2e4*0.85; 0; 0 ; 5e4 ; 0; 0.05; 1];
L = 0;
t0 = 3000;
t = 480;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
L = 1e3;
[T,X] = BennettModel(K,t,dt,L,x01);
% plot the time series of receptor species and their conservation
figure(4);
% plot time series of unphosphorylated surface receptors
subplot(1,3,1);
plot(T,X(:,1));
title('Unphosphorylated surface receptors');
xlabel('Time (seconds)');
ylabel('Number');
% plot time series of phosphorylated surfae receptors
subplot(1,3,2);
plot(T,X(:,2));
title('Phosphorylated surface receptors');
xlabel('Time (seconds)');
ylabel('Number');
% plot the summation of the two variables and its comparison to the total
% number of rceptors
subplot(1,3,3);
plot(T,X(:,1) + X(:,2),'b',T,ones(1:length(T),1).*K(1),'r');
title('Surface receptors');
xlabel('Time (seconds)');
ylabel('Number');
legend('R + Rp','Rtot');


%% plot G protein amount
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
%x0 = [2e4*0.85; 0; 0 ; 5e4 ; 0; 0.05; 1];
L = 0;
t0 = 3000;
t = 480;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
% Run simulation
L = 1e3;
[T,X] = BennettModel(K,t,dt,L,x01);
figure(5);
plot(T,X(:,3),'b',T,ones(1:length(T), 1).*K(8),'r');
title('G protein');
xlabel('Time (seconds)');
ylabel('Number');
legend('activated G protein','Gtot');

%% plot the time series of ip3, pip2, h and their conservation
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
%x0 = [2e4*0.85; 0; 0 ; 5e4 ; 0; 0.05; 1];
L = 0;
t0 = 3000;
t = 480;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
% Run simulation
L = 1e3;
[T,X] = BennettModel(K,t,dt,L,x01);
figure(6);
plot(T,X(:,4),'b',T,X(:,5)*K(18)*K(17),'g',T,ones(1:length(T), 1).*K(12),'r');
title('PIP2 associates');
xlabel('Time (seconds)');
ylabel('Number');
legend('PIP2','IP3','PIP2tot');
figure(7);

%% plot the time series of calcium species and their conservation



