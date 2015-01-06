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

        
%% init varaibles
GoodModels = cell(100,1); 

%%
% *For initial cell, perform ABC:*

%% determine parameters and functions
Eps=1.25;
Iter=50;
PopSize=24;

fRand = @(N) (rand(N,32)-0.5)*4;
fEval=@(Score,Pop) Score<Eps & max(abs(Pop),[],2)<4;
fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*BMparams,L,tInterp,CaFlt(:,1));

%% Find initial population


%% perform 50 iterations to sample from posterior
for i=1:100
    %%
    fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*BMparams,L,tInterp,CaFlt(:,i));
    FirstGenPop = rejectionSampling(fScore,fRand,fEval,PopSize);
    FirstGenWeights = ones(1,PopSize)/PopSize;
    GoodModels{i} = ABCrunGivenInitialPopulation(FirstGenPop,FirstGenWeights,fScore);
end

%% do for all models 
for i=1:100
    fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*BMparams,L,tInterp,CaFlt(:,i));
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