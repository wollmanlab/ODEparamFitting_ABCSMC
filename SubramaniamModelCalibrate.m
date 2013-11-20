%% This script tests the model with different combinations of parameters
% specify the parameters and ICs to be tweaked 
% icVec has the indices in the IC vector to be calibrated
icVec = [2];
%icRangeVec specify the ranges of the ICs
icRangeVec  = {4.11e-2};
% pVec has the indices in the parameter vector to be callibrated
pVec = [];
%pRangeVec specifies the ranges of the parameters to be calibrated
pRangeVec = {};
% Combine the two cell array of ranges into one cell array of ranges 
rangeVec = [icRangeVec pRangeVec];
%%
% make the multidimensional meshgrid
inputGrid = cell(1,length(icRangeVec) + length(pRangeVec));
[inputGrid{:}] = ndgrid(rangeVec{:});
% make the inputGrid into a matrix of input pairs
inputMatrix = zeros(size(inputGrid{:},1) + size(inputGrid{:},2),length(inputGrid));
for i = 1:size(inputMatrix,2)
    currInputGrid = inputGrid{:};
    inputMatrix(:,i) = currInputGrid(:);
end
%%
% go through the pairs and calculate hte fitness of cytosolic calcium
fitnessMatrix = zeros(size(inputMatrix,2),1);
for i = 1:size(inputMatrix,2)
   currInputVec = inputMatrix(i,:);
   ver = 0;
   [K,x0,L,t0,t,dt] = SubramaniamInput(ver);
   % change the ICs
   for col = 1:length(icVec)
       x0(icVec(col)) = currInputVec(col);
   end
   % chagne the parameters
   for col = 1:length(pVec)
       K(pVec(col)) = currInputVec(length(icVec) + col);
   end
   [T0L,X0L,speciesArray0,fluxMatrix0] = SubramaniamModel(K,t0,x0,dt);
   xL = X0L(end,:);
   xL(1) = L;
   [TL,XL,speciesArray1,fluxaMatrix1] = SubramaniamModel(K,t,xL,dt);
end


