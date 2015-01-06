clear
close all
clc
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
clf;
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
clf;
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
clf;
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
clf;
plot(T,X(:,4),'b',T,X(:,5)*K(18)*K(17),'g',T,ones(1:length(T), 1).*K(12),'r');
title('PIP2 associates');
xlabel('Time (seconds)');
ylabel('Number');
legend('PIP2','IP3','PIP2tot');

figure(7);
clf;
plot(T,X(:,7));
title('fraction of free IP3 receptor ');
xlabel('Time (seconds)');
ylabel('fraction');

%% plot the time series of responses and their variation with ligand cnocentration
lArray = [0; 1; 10; 1e2; 1e3 ];
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
tCell = cell(length(lArray),1);
xCell = cell( length(lArray),1);
for i = 1:length(lArray)
    L = 0;
    t0 = 3000;
    t = 240;
    % Run the simulation to equilibrium
    [T0,X0] = BennettModel(K,t0,dt,L,x0);
    x01 = X0(end,:)';
    % Run simulation
    L = lArray(i);
    [T,X] = BennettModel(K,t,dt,L,x01);
    tCell{i} = T; 
    xCell{i} = X;
end
clr = jet(5);
figure(8);
clf;
hold on;
title('Calcium Response Variation with varing ligand concentration');
xlabel('Time (seconds)');
ylabel('micromlar');
for i = 1:length(lArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,6),'color',clr(i,:));  
end
legend('0', '1', '10', '1e2', '1e3');
%% Plot the calcium responses to varying amounts PIP2 replenishment rate
rrArray = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10];
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
tCell = cell(length(rrArray),1);
xCell = cell( length(rrArray),1);
for i = 1:length(rrArray)
    L = 0;
    t0 = 3000;
    t = 240;
    % Run the simulation to equilibrium
    [T0,X0] = BennettModel(K,t0,dt,L,x0);
    x01 = X0(end,:)';
    % Run each of the simulation
    L = 1e3;
    K(13) = rrArray(i);
    [T,X] = BennettModel(K,t,dt,L,x01);
    tCell{i} = T; 
    xCell{i} = X;
end
clr = jet(length(rrArray));
figure(9);
clf;
hold on;
title('Calcium response to varying rr');
xlabel('Time (sec)');
ylabel('micromolar');
for i = 1:length(rrArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,6),'color',clr(i,:));  
end
legend(strread(num2str(rrArray),'%s'));
%% Plot the calcium responses to varying a2 parameter
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
a2Array = [K(24)*1e-3,K(24)*1e-2,K(24)*1e-1,K(24)*1, K(24)*10, K(24)*100];
tCell = cell(length(a2Array),1);
xCell = cell( length(a2Array),1);
for i = 1:length(a2Array)
    L = 0;
    t0 = 3000;
    t = 240;
    % Run the simulation to equilibrium
    [T0,X0] = BennettModel(K,t0,dt,L,x0);
    x01 = X0(end,:)';
    % Run each of the simulation
    L = 1e3;
    K(24) = a2Array(i);
    [T,X] = BennettModel(K,t,dt,L,x01);
    tCell{i} = T; 
    xCell{i} = X;
end
clr = jet(length(a2Array));
figure(10);
clf;
subplot(1,2,1);
hold on;
title('Calcium response to varying a2');
xlabel('Time (sec)');
ylabel('micromolar');
for i = 1:length(a2Array)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,6),'color',clr(i,:));  
end
legend(strread(num2str(a2Array),'%s'));
subplot(1,2,2);
hold on;
title('IP3 receptor fraction resposne to varying a2');
xlabel('Time (sec)');
ylabel('fraction');
for i = 1:length(a2Array)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,7),'color',clr(i,:));  
end
legend(strread(num2str(a2Array),'%s'));
%% Plot the calcium responses to varying kr parameter
ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
krArray = [K(4)*1e-3,K(4)*1e-2,K(4)*1e-1,2*K(4)*1e-1,3*K(4)*1e-1,4*K(4)*1e-1,5*K(4)*1e-1,6*K(4)*1e-1,7*K(4)*1e-1,8*K(4)*1e-1,9*K(4)*1e-1,K(4)*1,  K(4)*2,  K(4)*3, K(4)*4, K(4)*5,K(4)*10, K(4)*100];
tCell = cell(length(krArray),1);
xCell = cell( length(krArray),1);
for i = 1:length(krArray)
    L = 0;
    t0 = 3000;
    t = 1000;
    % Run the simulation to equilibrium
    [T0,X0] = BennettModel(K,t0,dt,L,x0);
    x01 = X0(end,:)';
    % Run each of the simulation
    L = 1e3;
    K(4) = krArray(i);
    [T,X] = BennettModel(K,t,dt,L,x01);
    tCell{i} = T; 
    xCell{i} = X;
end
clr = jet(length(krArray));
figure(11);
clf;
subplot(1,3,1);
hold on;
title('Calcium response to varying kr');
xlabel('Time (sec)');
ylabel('micromolar');
for i = 1:length(krArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,6),'color',clr(i,:));  
end
legend(strread(num2str(krArray),'%s'));
subplot(1,3,2);
hold on;
title('IP3 receptor fraction resposne to varying kr');
xlabel('Time (sec)');
ylabel('fraction');
for i = 1:length(krArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,7),'color',clr(i,:));  
end
legend(strread(num2str(krArray),'%s'));
subplot(1,3,3);
hold on;
title('unphosphorylated surface receptor');
xlabel('Time (sec)');
ylabel('fraction');
for i = 1:length(krArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,1),'color',clr(i,:));  
end
legend(strread(num2str(krArray),'%s'));

%% Plot the calcium responses to varying kp (phosphorylation rate)

ver = 1;
[K,t0,t,dt,x0] = BennettModelInput(ver);
kpArray = [K(5)*1e-3,K(5)*1e-2,K(5)*1e-1,2*K(5)*1e-1,3*K(5)*1e-1,4*K(5)*1e-1,5*K(5)*1e-1,6*K(5)*1e-1,7*K(5)*1e-1,8*K(5)*1e-1,9*K(5)*1e-1,K(5)*1, K(5)*10, K(5)*100];
tCell = cell(length(kpArray),1);
xCell = cell( length(kpArray),1);
for i = 1:length(kpArray)
    L = 0;
    t0 = 3000;
    t = 1000;
    % Run the simulation to equilibrium
    [T0,X0] = BennettModel(K,t0,dt,L,x0);
    x01 = X0(end,:)';
    % Run each of the simulation
    L = 1e3;
    K(5) = kpArray(i);
    [T,X] = BennettModel(K,t,dt,L,x01);
    tCell{i} = T; 
    xCell{i} = X;
end
clr = jet(length(kpArray));
figure(12);
clf;
subplot(1,3,1);
hold on;
title('Calcium response to varying kp');
xlabel('Time (sec)');
ylabel('micromolar');
for i = 1:length(kpArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,6),'color',clr(i,:));  
end
legend(strread(num2str(kpArray),'%s'));
subplot(1,3,2);
hold on;
title('IP3 receptor fraction resposne to varying kp');
xlabel('Time (sec)');
ylabel('fraction');
for i = 1:length(kpArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,7),'color',clr(i,:));  
end
legend(strread(num2str(kpArray),'%s'));
subplot(1,3,3);
hold on;
title('unphosphorylated surface receptor');
xlabel('Time (sec)');
ylabel('fraction');
for i = 1:length(kpArray)
  currXCell = xCell{i};
  plot(tCell{i},currXCell(:,1),'color',clr(i,:));  
end
legend(strread(num2str(kpArray),'%s'));
%% 
ver = 3;
[K,t0,t,dt,x0] = BennettModelInput(ver);
K(4) =0;
L = 0;
t0 = 3000;
t = 900;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
% Run each of the simulation
LArray = [1e2, 1e1,8, 6, 4, 2 ,1,0.8, 0.4, 0.1,0.01];
xArray = cell(length(LArray),1); 
tArray = cell(length(LArray),1);
for i =1:length(LArray)
    L = LArray(i);
    [tArray{i},xArray{i}] = BennettModel(K,t,dt,L,x01);
end
clr = jet(length(LArray));
figure(1);
clf;
subplot(1,3,1);
hold on
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,6),'color',clr(i,:));
end
title('calcium concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))),strcat('L = ',num2str(LArray(7))), strcat('L = ',num2str(LArray(8))), strcat('L = ',num2str(LArray(9))),strcat('L = ',num2str(LArray(10))), strcat('L = ',num2str(LArray(11))) );
subplot(1,3,2);
hold on
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,2),'color',clr(i,:));
end
title('phosphorylated surface receptor with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))),strcat('L = ',num2str(LArray(7))), strcat('L = ',num2str(LArray(8))), strcat('L = ',num2str(LArray(9))),strcat('L = ',num2str(LArray(10))), strcat('L = ',num2str(LArray(11))) );
subplot(1,3,3);
hold on
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,1),'color',clr(i,:));
end
title('phosphorylated surface receptor with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))),strcat('L = ',num2str(LArray(7))), strcat('L = ',num2str(LArray(8))), strcat('L = ',num2str(LArray(9))),strcat('L = ',num2str(LArray(10))), strcat('L = ',num2str(LArray(11))) );
%% Plot the estimate of LRs/(K1 + L)
ver = 3;
[K,t0,t,dt,x0] = BennettModelInput(ver);
L = 0;
t0 = 3000;
t = 900;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
% Run each of the simulation
LArray = [1e2, 1e1, 1,0.1,0.01];
xArray = cell(length(LArray),1); 
tArray = cell(length(LArray),1);
for i =1:length(LArray)
    L = LArray(i);
    [tArray{i},xArray{i}] = BennettModel(K,t,dt,L,x01);
end
clr = jet(length(LArray));
figure(2);
clf;
subplot(1,2,1);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},LArray(i)*X(:,1)./(K(2)+ LArray(i)),'color',clr(i,:));
end
title('LRs/(k1 + L)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))) );
% get the value of kp*L/(K1 + L)Rs
subplot(1,2,2);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},K(5)*LArray(i)*X(:,1)./(K(2)+ LArray(i)),'color',clr(i,:));
end
title('KpLRs/(k1 + L)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))) );

%% Plot the modification of Bennett model with varying ligand concentration
ver = 3;
[K,t0,t,dt,x0] = BennettModelInput(ver);
L = 0;
t0 = 3000;
t = 1000;
kpp = 50;
km = 400;
%K(16) = K(16)*100;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t0,dt,L,x0);
x01 = X0(end,:)';
% Run each of the simulation
%LArray = [1e2, 1e1,8,6,4, 1,0.8, 0.6, 0.4,0.1,0.01];
LArray = [10,3.33,1,0.33,0.1,0 ];
xArray = cell(length(LArray),1); 
tArray = cell(length(LArray),1);

for i =1:length(LArray)
    L = LArray(i);
    [tArray{i},xArray{i}] = BennettModelMod(K,t,dt,L,x01,kpp,km);
end
clr = jet(length(LArray));
figure(1);
clf;
hold on
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,6),'color',clr(i,:));
end
title('calcium concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );


figure(3);
clf;
subplot(1,3,1);
hold on
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,6),'color',clr(i,:));
end
title('calcium concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );

subplot(1,3,2);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,1),'color',clr(i,:));
end
title('R unphosphorylated concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))),strcat('L = ',num2str(LArray(6))) );

subplot(1,3,3);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,2),'color',clr(i,:));
end
title('R phosphorylated concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );

figure(4);
subplot(2,2,1);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,3),'color',clr(i,:));
end
title('Activated G protein concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );

subplot(2,2,2);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,4),'color',clr(i,:));
end
title('PIP2 concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))) , strcat('L = ',num2str(LArray(6))) );

subplot(2,2,3);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,5),'color',clr(i,:));
end
title('IP3 concentration with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))) , strcat('L = ',num2str(LArray(6))) );

subplot(2,2,4);
hold on;
for i = 1:length(LArray)
    X = xArray{i};
    plot(tArray{i},X(:,7),'color',clr(i,:));
end
title('Active IP3R with varying ligand');
ylabel('concentration(micromolar)');
xlabel('time (seconds)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );
%% Try to vary the set of parameters over a certain range to see if it will affect the calcium profile
range = [5 6 7 8 9 10];
% indicate the parameter to be tweaked
paramInd = 9;
% Establish the list of ligand values as perturbation
LArray = [10,3.33,1,0.33,0.1,0 ];
% Determine the number of subplots that need to be made
m =2;
n = ceil(length(range)/2);
% go through the list of parameter tweak vallues
figure;
clf;
for i = 1:length(range)
    % Get the new inputs
    % The parameter tweak of the current modification of the model
    ver = 3;
    [K,t0,t,dt,x0] = BennettModelInput(ver);
    % change the parameter according to the paramList
    paramList = range.*K(paramInd);
    K(2) = 10;
    K(paramInd) = paramList(i);
    L = 0;
    t0 = 3000;
    t = 1000;
    kpp = 50;
    km = 400;
    % Run the simulation to equilibrium
    %[T0,X0] = BennettModel(K,t0,dt,L,x0);
    [T0,X0] = BennettModelMod(K,t0,dt,L,x0,kpp,km);
    x01 = X0(end,:)';
    subplot(m,n,i);
    hold on;
    xArray = cell(length(LArray),1); 
    tArray = cell(length(LArray),1);
    % Go through the elements of LArray to form one subplot
    for j =1:length(LArray)
        L = LArray(j);
        [tArray{j},xArray{j}] = BennettModelMod(K,t,dt,L,x01,kpp,km);
    end
    clr = jet(length(LArray));
    % plot the current subplot
    for ind = 1:length(LArray)
        X = xArray{ind};
        plot(tArray{ind},X(:,6),'color',clr(ind,:));
    end
    title(strcat('Parameter K(',num2str(paramInd),' ) = ',num2str(paramList(i))))
    ylabel('concentration(micromolar)');
    xlabel('time (seconds)');
    legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );
end
%% plot the experimental plot and the computational model plot together
ver = 3;
[K,t0,t,dt,x0] = BennettModelInput(ver);
K(2) =10;
K(9) = 11.25;
L = 0;
t0 = 3000;
t = 1000;
kpp = 50;
km = 400;
% Run the simulation to equilibrium
[T0,X0] = BennettModelMod(K,t0,dt,L,x0,kpp,km);
x01 = X0(end,:)';
LArray = [10,3.33,1,0.33,0.1,0 ];
xArray = cell(length(LArray),1); 
tArray = cell(length(LArray),1);
for j =1:length(LArray)
    L = LArray(j);
    [tArray{j},xArray{j}] = BennettModelMod(K,t,dt,L,x01,kpp,km);
end
clr = jet(length(LArray));
figure(2)
%subplot(1,2,2)
%clf;
hold on;
% plot the current subplot
for ind = 1:length(LArray)
    X = xArray{ind};
    time = tArray{ind};
    %plot(time(1:500),X(1:500,6),'color',clr(ind,:));
    plot(time(1:500),X(1:500,6)/X(1,6),'color',clr(ind,:),'LineWidth',6);
end
title('Computational Modeling')
ylabel('Calcium Fold Increase');
xlabel('Time (sec)');
legend(strcat('L = ',num2str(LArray(1))), strcat('L = ',num2str(LArray(2))), strcat('L = ',num2str(LArray(3))), strcat('L = ',num2str(LArray(4))), strcat('L = ',num2str(LArray(5))), strcat('L = ',num2str(LArray(6))) );
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',20,'fontWeight','bold');
set(gca,'fontsize',20);
%% 3/17/2014
%{
   Objective: test the hypothesis that the cell-to-cell variation in the PLC flux
   related parameters is sufficient to cause the bnary calcium response in
   the cells perturbed by ATP
%}
% generate input
close all; clear all;
ver = 3;
[K,t0,t,dt,x0] = BennettModelInput(ver);
K(2) =10;
K(9) = 11.25;
t0 = 3000;
t = 1000;
kpp = 50;
km = 400;
paramRange = [0.1, 0.33,0.66,1,3.33,6.66,10];
% vary alpha, the effectve sgnal gain parameter in the equation for rh
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(16) = KMod(16)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(1);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying alpha');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end

% vary Kc, the dissociation constant for Ca binding site on the PLC
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(15) = KMod(15)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(2);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying Kc');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end




% vary xi, the fraction of mobile receptors
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(7) = KMod(7)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(3);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying xi');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        %legend('0.1','0.5','1','5','10');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end





% vary K1, the 
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(2) = KMod(2)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(4);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');       
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
       legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying K1');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end





% vary delta, the ratio of activities of the ligand unbound and bound
% receptor species
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(14) = KMod(14)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(5);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');  
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying delta');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end









% vary ka, the G protein activation rate parameters
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(10) = KMod(10)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(6);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');       
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying ka');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end










% vary kd, the G protein deactivation rate parameters
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(11) = KMod(11)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(7);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');     
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying kd');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end




% vary Rt, the total number of receptor
xArray = cell(length(paramRange),1);
tArray = cell(length(paramRange),1);
x0Array = cell(length(paramRange),1);
t0Array = cell(length(paramRange),1);
for j =1:length(paramRange)
    % Run the simulation to equilibrium
    KMod = K;
    KMod(1) = KMod(1)*paramRange(j);
    L = 0;
    [t0Array{j},x0Array{j}] = BennettModelMod(KMod,t0,dt,L,x0,kpp,km);
    X0 = x0Array{j};
    x01 = X0(end,:)';
    % set the ligand concentration, the range of the parameter, and the indices
    % of parameters
    L = 10;
    [tArray{j},xArray{j}] = BennettModelMod(KMod,t,dt,L,x01,kpp,km);
end
figure(8);
clf;
hold on;
clr = jet(length(paramRange));
% plot all state variables
for i = 1:7
    % unphosphorylated surface receptor
    if i ==1
        subplot(2,4,i);
        hold on; 
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            %x0 = x0Array{ind};
            %t0 = t0Array{ind};
            %plot(t0,x0(:,i),'color',clr(ind,:))
            %plot(max(t0)+time,X(:,i),'color',clr(ind,:));
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Unphos. R, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % phosphorylated surface receptor
    elseif i ==2
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('Phos. R, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');      
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');       
    % activated G protein    
    elseif i ==3
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('G protein, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  PIP2
    elseif i ==4
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('PIP2, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    %  IP3   
    elseif i ==5
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('IP3, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % cytosolic calcuim     
    elseif i ==6
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('calcium, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    % h, fraction of activated IP3R     
    elseif i ==7
        subplot(2,4,i);
        hold on;
        for ind = 1:length(paramRange)
            X = xArray{ind};
            time = tArray{ind};
            plot(time,X(:,i),'color',clr(ind,:));
        end
        title('h, varying Rt');
        ylabel('concentration(micromolar)');
        xlabel('time (seconds)');        
        legend('0.1', '0.33','0.66','1','3.33','6.66','10');
    end
end



