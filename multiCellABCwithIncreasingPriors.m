%% Init
clear 
close all
clc

warning off

%% include
addpath /home/jason/calciummodeling/
addpath /home/jason/calciummodeling/Traj_Compare_Conditions/Traj_Orig_Ca/

%% load experimental data
Ca=[]; 
CellId={}; 
CaFlt={}; 
L=[]; 
load HundredRandomCells 

%% determine the "base" parameters 
BMparams=BennettModelInput(1); 
BMparams(17:19)=[]; 
BMparams=BMparams(:)';

RoyParamSet=BMparams.*10.^[1.2096 -0.7784 -3.4222  2.2206 -0.3126 -0.0611  1.0283 1.1984 ...
             0.9949 -2.6569 -2.1135 -0.8790 -0.5521 -1.9549  2.7077 1.6535 ...
             1.1794  0.4340  0.6620  1.7681 -1.3118  0.9545 -1.5936 1.5925 ...
            -0.4810 -0.2436 -2.6716 -0.0377 -0.1247  0.6700 -1.2272 0.7170]; 

        
%% init varaibles
GoodModels = cell(100,1); 

%%
% *For initial cell, perform ABC:*

%% decide which cells is the "most central"
D=squareform(pdist(CaFlt'));
[~,cellid]=min(sum(D)); 

%% determine parameters and functions
Eps=1.25;
Iter=50;
PopSize=24;

fRand = @(N) (rand(N,32)-0.5)*4;
fEval=@(Score,Pop) Score<Eps & max(abs(Pop),[],2)<4;
fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*RoyParamSet,L,tInterp,CaFlt(:,1));

%% Find initial population
FirstGenPop = rejectionSampling(fScore,fRand,fEval,PopSize);
FirstGenWeights = ones(1,PopSize)/PopSize;

%% perform 50 iterations to sample from posterior
for i=1:100
    fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*RoyParamSet,L,tInterp,CaFlt(:,i));
    GoodModels{i} = ABCrunGivenInitialPopulation(FirstGenPop,FirstGenWeights,fScore);
end

%% do for all models 
for i=1:100
    fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*RoyParamSet,L,tInterp,CaFlt(:,i));
    for j=1:size(GoodModels{i},1)
        [F{i}(j),Sim{i}{j},F0{i}(j),F1{i}(j)]=fScore(GoodModels{i}(j,:)); 
    end
    i
end

% %% Continue for all cells
% for cellid=2:size(CaFlt,2)
%     
%     %% update the scoring function based on new cell data
%     fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*RoyParamSet,L,tInterp,CaFlt(:,cellid));
%     
%     %% determine initial sample based on weights with within covariance as initial sample. 
%     AllModels = cat(1,GoodModels{:}); 
%     WithinCov = cellfun(@cov,GoodModels(1:cellid-1),'uniformoutput',0); 
%     WithinCov = mean(cat(3,WithinCov{:})); 
%     
%     CrossScores = fScore(AllModels);
%     
%     fRand = @(N) AllModels(randsample(AllModels,N,true,CrossScores),:) + ...
%                  mvnrnd(zeros(1,32),2*WithinCov,N);
%              
%     %% Find initial population for cell cellid
%     FirstGenPop = rejectionSampling(fScore,fRand,fEval,PopSize);
%     FirstGenWeights = ones(1,PopSize)/PopSize; 
%     
%     %% find good models
%     GoodModels{cellid} = ABCrunGivenInitialPopulation(FirstGenPop,FirstGenWeights,fScore);
% 
%     %% save results
%     save HundredCellFitting
%     
% end
% 
%     
% %% save results
% save HundredCellFitting