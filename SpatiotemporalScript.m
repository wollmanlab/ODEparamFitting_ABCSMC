%{
   This is a script that will eventually be turned into a function.
   It simulates the spatiotemporal behavior of ATP for right now. It will
   eventually be incorporated with the intracellular model of calcium. 
   
   Input:
         D = diffusion coefficient of ATP
         s = diameter of the wound size
         numATP = total number of ATP storage per cell
         h = height of cylinder
         numCell = number of cells in the simulation
         xDim = x dimension of the simulation 
         yDim = y dimension of the simulation 
         dt = size of time step
         dDim = unit of spatial increment
         mode = mode of ATP release. 1 = dilution, 2 = diffusion with
         degradation, 3 = wave propagation with positive feedback
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
%}
%% First make the grid space of the simulaiton
% generate inputs
ver = 0;
[s,height,D,numATP,numCell,r,dr,tNum,dt,mode] = SpatioTemporalInput(ver);
IC = zeros(r,1);
IC(1:2) = numATP; % for current testing purpose just use this IC
%% Check for the mode of propagation 
switch mode
    % dilution
    case 1
        m = RadialDiffusion(D,IC,dr,dt,tNum,1);
        
    % diffusion with degradation
    case 2
        
    % wave propagation
    case 3
        
end
%% plot the results in ten snapshots
clr = jet(11);
figure(1);
clf;
hold on;
rVec = dr:dr:r*dr;
%plot(rVec,m(:,1),'color',clr(1,:));
for i = 1:10
    plot(rVec,m(:,i*tNum/10+1),'color',clr(i+1,:));
end
title('Ligand Diffusion');
xlabel('Radial distance (micrometer)');
ylabel('# Molecules');
legend(strcat('t = ',num2str((tNum/10)*dt)), strcat('t = ',num2str((2*tNum/10)*dt)), strcat('t = ',num2str((3*tNum/10)*dt)), strcat('t = ',num2str((4*tNum/10)*dt)), strcat('t = ',num2str((5*tNum/10)*dt)), strcat('t = ',num2str((6*tNum/10)*dt)), strcat('t = ',num2str((7*tNum/10)*dt)), strcat('t = ',num2str((8*tNum/10)*dt)), strcat('t = ',num2str((9*tNum/10)*dt)), strcat('t = ',num2str((10*tNum/10)*dt))); 
%% plot the conservation
figure(2);
mat = m(2:end,:); 
ligandSum = sum(mat,1);
plot(0:dt:tNum*dt, ligandSum);
title('ligand conservation');
xlabel('Time(sec)');
ylabel('Sum of ligand');