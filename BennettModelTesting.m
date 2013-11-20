%% This script tests the accuracy of the Bennett Model by reproducing the plots in the paper
%% Reproduce Fig.3 (a) in the paper that relates ligand concentration to the percentage of surface receptors at equilibrium
% 
[K,t0,t,dt,x0] = BennetModelInput();
L = 0;
t0 = 3000;
t = 2100;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t,dt,L,x0);
x0 = X0(end,:)';
receptorArray = zeros(5,1); 
utpArray = [10^-1; 1; 10; 10^2 ; 10^3];
% Get data points for the receptor array to be plotted
for i = 1:size(receptorArray,1)
    [T,X] = BennettModel(K,t,dt,utpArray(i),x0);
    receptorArray(i) = (X(end,1) + X(end,2))/K(1);
end
figure(1);
clf;
subplot(1,2,1);
plot(utpArray,receptorArray);
title('Surface Receptors at Equilibrium');
xlabel('ligand (micromolar)');
ylabel('Surface Receptors Percentage');


L = 10^3;
[T,X] = BennettModel(K,t,dt,L,x0);
subplot(1,2,2);
surfRecep = (X(:,1) + X(:,2))/K(1);
plot(T,surfRecep);
title('Surface receptor time course');
xlabel('Time (second)');
ylabel('Surface receptor( fraction)');


%% Reproduce Fig. 5 in the Bennett paper
[K,t0,t,dt,x0] = BennetModelInput();
L = 0;
t0 = 3000;
t = 480;
% Run the simulation to equilibrium
[T0,X0] = BennettModel(K,t,dt,L,x0);
x0 = X0(end,:)';
rrArray = [10;0.1;0.015];
tArray = cell(3,1);
xArray = cell(3,1);
% Perform simulations
for i = 1:size(rrArray,1)
    K1 = K;
    K1(13) = rrArray(i);
    L = 10^3;
    [T,X] = BennettModel(K,t,dt,L,x0);
    tArray{i} = T; 
    xArray{i} = X;
end
% Plot the figures
figure(2);
% plot Figure 5 (a)
subplot(2,2,1);
xArray1 = xArray{1};
xArray2 = xArray{2};
xArray3 = xArray{3};
plot(tArray{1}/2,xArray1(:,5),'r',tArray{2}/2,xArray2(:,5),'b',tArray{3}/2,xArray3(:,5),'g');
title( 'IP3 variaiton with rr');
xlabel('Time (seconds)');
ylabel('IP3 (micromolar)');
legend('rr = 10','rr = 0.1','rr = 0.015');
% plot Figure 5(b)
subplot(2,2,2);
plot(tArray{1}/2,xArray1(:,6),'r',tArray{2}/2,xArray2(:,6),'b',tArray{3}/2,xArray3(:,6),'g');
title(' Calcium variation with rr');
xlabel('Time (seconds)');
ylabel('Calcium (micromolar)');
legend('rr = 10','rr = 0.1', 'rr = 0.015');
% plot figure 5(c)
subplot(2,2,3);
plot(tArray{1},xArray1(:,4)/K(12),'r',tArray{2},xArray2(:,4)/K(12),'b',tArray{3},xArray3(:,4)/K(12),'g');
title(' PIP2 variation with rr');
xlabel('Time (seconds)');
ylabel('PIP2 (micromolar)');
legend('rr = 10','rr = 0.1', 'rr = 0.015');
% plot figure 5 (d)
subplot(2,2,4);
plot(tArray{1}/2,xArray1(:,4)/K(8));
title('Fraction of activated G protein');
xlabel('Time (seconds)');
ylabel('Fraction (micromolar)');
%%
