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
%%
clear
close all
clc
stkshow('close','all')
%% First make the grid space of the simulaiton
% generate inputs
[dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(300,3000);

%% Use the analytical solution to draw a comparison to to make sure there is no artifacts introduced by numerics
%plot the comparison between the finite difference implementation and the analytical expresison
% define the 2 dimensional 
tVec = 0:dt:tNum*dt;
if mod(max(xDim,yDim),2) == 0
    rDim = max(xDim,yDim)/2 +1;
else
   rDim = ceil(max(xdim,yDim)/2);
end
rVec = 0:dr:(rDim-1)*dr;
h = 100;
%funcC2d = @(r,t) totalLig*(dr^2)/(4*pi*D.*t)*exp(-r.^2/4/D./t)*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
%funcC2d = @(r,t) totalLig/(4*pi*D.*t)*exp(-r.^2/4/D./t)*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);;
funcC2d = @(r,t) totalLig./(4*pi.*D.*t./h).*exp(-r.^2/4/D./t)/0.6/1000;
figure(1);
clf;
subplot(1,2,1);
hold on;
clr = jet(10);
for i = 1:10
    plot(rVec,funcC2d(rVec,tVec(i*tNum/10 +1)),'color',clr(i,:));
end
title('Analytical solution of dilution');
xlabel('Radial distance(micrometer)');
ylabel('micromolar');
legend(strcat('t = ',num2str(tVec(tNum/10 +1))), strcat('t = ',num2str(tVec(2*tNum/10 +1))), strcat('t = ',num2str(tVec(3*tNum/10 +1))), strcat('t = ',num2str(tVec(4*tNum/10 +1))), strcat('t = ',num2str(tVec(5*tNum/10 +1))), strcat('t = ',num2str(tVec(6*tNum/10 +1))), strcat('t = ',num2str(tVec(7*tNum/10 +1))), strcat('t = ',num2str(tVec(8*tNum/10 +1))), strcat('t = ',num2str(tVec(9*tNum/10 +1))), strcat('t = ',num2str(tVec(10*tNum/10 +1)))  ); 
%%
% numerical simualtion
m1 = TwoDimDiffusion2(D,totalLig,xDim,yDim,dr,dt,tNum,1,deg,reaction);
clr = jet(10);
subplot(1,2,2);
hold on;
increment = tNum/10;
conVec = zeros(1,1);
%%
for i = 1:10
        % get the current snap shot
        currentM = m1{i*increment+1};
        % The center of diffusion
        centerP = [ceil(size(currentM,1)/2),ceil(size(currentM,2)/2)];
        % the vector in the radial direction from the center
        conVec = currentM(centerP(1),centerP(2):size(currentM,2));
        % convert conVec to units of concentration in micromolar(divided by 100 micron height)
        conVec = conVec*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
        % radial position vector
        rVec = 0:dr:dr*(ceil(size(currentM,1)/2)-1);
        plot( rVec,conVec,'color',clr(i,:));
end
%%
for i = 1:10
        % get the current snap shot
        conVec = m1{i*increment+1};
        % convert conVec to units of concentration in micromolar(divided by 100 micron height)
        conVec = conVec*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
        % radial position vector
        rVec = 0:dr:dr*(ceil(size(currentM,1)/2)-1);
        plot( rVec,conVec,'color',clr(i,:));
end

%%
m2 = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,2,deg,reaction);
clr = jet(10);
figure(3);
clf;
hold on;
increment = (tNum)/10;
for i = 1:10
        % get the current snap shot
        currentM = m2{i*increment+1};
        % The center of diffusion
        centerP = [ceil(size(currentM,1)/2),ceil(size(currentM,2)/2)];
        % the vector in the radial direction from the center
        conVec = currentM(centerP(1),centerP(2):size(currentM,2));
        % convert conVec to units of concentration in micromolar(divided by 100 micron height)
        conVec = conVec*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
        % radial position vector
        rVec = 0:dr:dr*(ceil(size(currentM,1)/2)-1);
        plot( rVec,conVec,'color',clr(i,:));
end
%%
m3 = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,3,deg,reaction);
clr = jet(10);
figure(4);
clf;
hold on;
increment = (tNum)/10;
for i = 1:10
        % get the current snap shot
        currentM = m3{i*increment+1};
        % The center of diffusion
        centerP = [ceil(size(currentM,1)/2),ceil(size(currentM,2)/2)];
        % the vector in the radial direction from the center
        conVec = currentM(centerP(1),centerP(2):size(currentM,2));
        % convert conVec to units of concentration in micromolar(divided by 100 micron height)
        conVec = conVec*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
        % radial position vector
        rVec = 0:dr:dr*(ceil(size(currentM,1)/2)-1);
        plot( rVec,conVec,'color',clr(i,:));
end

%}
%% plan the series of computational experiments
% hypothesis 1 with varying wound sizes
hypo1WArray = cell(7,1);
% hypothesis 2 with varying wound sizes 
hypo2WArray = cell(7,1);
% hypothesis 3 with varying wuond sizes
hypo3WArray = cell(7,1);
% hypothesis 1 with varying cell density
hypo1DArray = cell(9,1);
% hypothesis 2 with varying cell density
hypo2DArray = cell(9,1);
% hypothesis 3 with varying cell density
hypo3DArray = cell(9,1);

%% Simulate hypotheses with varying wound sizes
woundSizeArray = [150, 200,250, 300,350,400, 450 ];
for i = 1:length(hypo1WArray)
    % generate the input sets
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(woundSizeArray(i),3000); 
    m = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,1);
    % Calculate the radial concentration
    centerP = [ceil(size(m{1},1)/2),ceil(size(m{1},2)/2)];
    radialCon = zeros(length(m),size(m{1},2)-centerP(2)+1);
    % loop through the cell array m to generate concentraiton matrix
    for  ind = 1:length(m) 
        % the current concentraiton matrix
        currM = m{ind};
        % read the radial concentration and convert to 
        radialCon(ind,:)  = currM(centerP(1),centerP(2):size(currM,2))*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
    end
    hypo1WArray{i} = radialCon;
end
% the threshold concentration for cell response to ligand. The units are
% ligand # per micron squared.
cellThres = 0.33; 
% construct the time vector
tVec = 0:dt:tNum*dt;
% loop through all the versions of inputs
% calculate the response distances
figure(2);
clf;
subplot(1,3,1);
hold on;
clr = jet(length(hypo1WArray));
% loop through the levels of wound sizes
for i = 1:length(hypo1WArray)
    % loop through the individual element radialConc matrices
    currentStack = hypo1WArray{i};
    % the radial Dist is the vector that describes the boundary that meets
    % the threshold concentration
    radialDist = zeros(length(tVec),1);
    % go through the time steps
    for ind  = 1:size(currentStack,1)
        radialPos = find(currentStack(ind,:)>= cellThres,1,'last');
        if isempty(radialPos)
           radialDist(ind) = 0; 
        else
           radialDist(ind) = radialPos;
        end
    end
    % interpolate the vector
    xq = 0:dt/10:tNum*dt;
    radialInterp = interp1(tVec,radialDist*dr,xq);
    plot(xq,radialInterp,'color',clr(i,:));
end
title('Threshold Distances under Dilution');
xlabel('Time (sec)');
ylabel('Threshold Distance');
legend('150','200','250','300','350','400','450');

% simulate hypothesis 2

for i = 1:length(hypo2WArray)
    % generate the input sets
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(woundSizeArray(i),3000); 
    m = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,2,deg,reaction);
    % Calculate the radial concentration
    centerP = [ceil(size(m{1},1)/2),ceil(size(m{1},2)/2)];
    radialCon = zeros(length(m),size(m{1},2)-centerP(2)+1);
    % loop through the cell array m to generate concentraiton matrix
    for  ind = 1:length(m) 
        % the current concentraiton matrix
        currM = m{ind};
        % read the radial concentration and convert to 
        radialCon(ind,:)  = currM(centerP(1),centerP(2):size(currM,2))*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
    end
    hypo2WArray{i} = radialCon;
end
% the threshold concentration for cell response to ligand. The units are
% ligand # per micron squared.
cellThres = 0.33; 
% construct the time vector
tVec = 0:dt:tNum*dt;
% loop through all the versions of inputs
% calculate the response distances
subplot(1,3,2);
hold on;
clr = jet(length(hypo2WArray));
% loop through the levels of wound sizes
for i = 1:length(hypo2WArray)
    % loop through the individual element radialConc matrices
    currentStack = hypo2WArray{i};
    % the radial Dist is the vector that describes the boundary that meets
    % the threshold concentration
    radialDist = zeros(length(tVec),1);
    % go through the time steps
    for ind  = 1:size(currentStack,1)
        radialPos = find(currentStack(ind,:)>= cellThres,1,'last');
        if isempty(radialPos)
           radialDist(ind) = 0; 
        else
           radialDist(ind) = radialPos;
        end
    end
    % interpolate the vector
    xq = 0:dt/10:tNum*dt;
    radialInterp = interp1(tVec,radialDist*dr,xq);
    plot(xq,radialInterp,'color',clr(i,:));
end
title('Threshold Distances under Degradation');
xlabel('Time (sec)');
ylabel('Threshold Distance');
legend('150','200','250','300','350','400','450');


%{
% simualate hypothesis 3
for i = 1:length(woundSizeArray)
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(woundSizeArray(i),3000);
    hypo2WArray{i} = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,3,deg,reaction);
end
%}
%% Simulate hypothesis with varying cell density
cellDensArray = [1/3,2/3, 1,4/3, 5/3, 2,7/3,8/3, 3 ]*3000;
% simulate hypothesis 1
for i = 1:length(hypo1DArray)
    % generate the input sets
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(300,cellDensArray(i)); 
    m = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,1);
    % Calculate the radial concentration
    centerP = [ceil(size(m{1},1)/2),ceil(size(m{1},2)/2)];
    radialCon = zeros(length(m),size(m{1},2)-centerP(2)+1);
    % loop through the cell array m to generate concentraiton matrix
    for  ind = 1:length(m) 
        % the current concentraiton matrix
        currM = m{ind};
        % read the radial concentration and convert to 
        radialCon(ind,:)  = currM(centerP(1),centerP(2):size(currM,2))*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
    end
    hypo1DArray{i} = radialCon;
end
% the threshold concentration for cell response to ligand. The units are
% ligand # per micron squared.
cellThres = 0.33; 
% construct the time vector
tVec = 0:dt:tNum*dt;
% loop through all the versions of inputs
% calculate the response distances
figure(3);
clf;
subplot(1,3,1);
hold on;
clr = jet(length(hypo1DArray));
% loop through the levels of wound sizes
for i = 1:length(hypo1DArray)
    % loop through the individual element radialConc matrices
    currentStack = hypo1DArray{i};
    % the radial Dist is the vector that describes the boundary that meets
    % the threshold concentration
    radialDist = zeros(length(tVec),1);
    % go through the time steps
    for ind  = 1:size(currentStack,1)
        radialPos = find(currentStack(ind,:)>= cellThres,1,'last');
        if isempty(radialPos)
           radialDist(ind) = 0; 
        else
           radialDist(ind) = radialPos;
        end
    end
    % interpolate the vector
    xq = 0:dt/10:tNum*dt;
    radialInterp = interp1(tVec,radialDist*dr,xq);
    plot(xq,radialInterp,'color',clr(i,:));
end
title('Threshold Distances under Dilution');
xlabel('Time (sec)');
ylabel('Threshold Distance');
legend('1/3','2/3','1','4/3','5/3','2','7/3','8/3','3');



% simulate hypothesis 2

for i = 1:length(hypo2DArray)
    % generate the input sets
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(300,cellDensArray(i)); 
    m = TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,2,deg,reaction);
    % Calculate the radial concentration
    centerP = [ceil(size(m{1},1)/2),ceil(size(m{1},2)/2)];
    radialCon = zeros(length(m),size(m{1},2)-centerP(2)+1);
    % loop through the cell array m to generate concentraiton matrix
    for  ind = 1:length(m) 
        % the current concentraiton matrix
        currM = m{ind};
        % read the radial concentration and convert to 
        radialCon(ind,:)  = currM(centerP(1),centerP(2):size(currM,2))*(1/(dr^2))*(1/100)*(1/(6.02e23))*(1e6)*(1e18)*(1/1e3);
    end
    hypo2DArray{i} = radialCon;
end
% the threshold concentration for cell response to ligand. The units are
% ligand # per micron squared.
cellThres = 0.33; 
% construct the time vector
tVec = 0:dt:tNum*dt;
% loop through all the versions of inputs
% calculate the response distances
subplot(1,3,2);
hold on;
clr = jet(length(hypo2DArray));
% loop through the levels of wound sizes
for i = 1:length(hypo2DArray)
    % loop through the individual element radialConc matrices
    currentStack = hypo2DArray{i};
    % the radial Dist is the vector that describes the boundary that meets
    % the threshold concentration
    radialDist = zeros(length(tVec),1);
    % go through the time steps
    for ind  = 1:size(currentStack,1)
        radialPos = find(currentStack(ind,:)>= cellThres,1,'last');
        if isempty(radialPos)
           radialDist(ind) = 0; 
        else
           radialDist(ind) = radialPos;
        end
    end
    % interpolate the vector
    xq = 0:dt/10:tNum*dt;
    radialInterp = interp1(tVec,radialDist*dr,xq);
    plot(xq,radialInterp,'color',clr(i,:));
end
title('Threshold Distances under degradation');
xlabel('Time (sec)');
ylabel('Threshold Distance');
legend('1/3','2/3','1','4/3','5/3','2','7/3','8/3','3');

%{
% simulate hypothesis 2
for i =1:length(cellDensArray)
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(300,cellDensArray(i));
    hypo1DArray{i} =  TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,2,deg);
end
% simulate hypothesis 3
for i =1:length(cellDensArray)
    [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(300,cellDensArray(i));
    hypo1DArray{i} =  TwoDimDiffusion(D,totalLig,xDim,yDim,dr,dt,tNum,3,deg,reaction);
end
%}


%% plot the heatmaps. Take wound size of 300 micron and density of 3000 cells per image
% plot heatmap of hypothesis 1 
figure(5);
clf;
% Get the radial concentration matrix for hypothesis 1, 300 micron wound
% size and 3000 per image
matrix1 = hypo1WArray{4};
clim = [0 1000];
tVec = 0:dt:tNum*dt;
rVec = 0:dr:dr*(size(matrix1,2)-1);
imagesc(rVec,tVec,matrix1,clim);
colormap(jet);
% plot heatmap of hypothesis 2

% plot heatmap of hypothesis 3

%% plot the maximum response distances vs. density and maximum response distances vs wound sizes
figure(4);

