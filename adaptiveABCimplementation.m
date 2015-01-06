%% Init
clear 
close all
clc

warning off

%% include
addpath /home/jason/calciummodeling/
addpath /home/jason/calciummodeling/Traj_Compare_Conditions/Traj_Orig_Ca/

%% load experimental data
load /home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/traj_147.mat
load /home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/tInterp.mat
L=10; 

% traj_interp_raw=traj_interp; 
traj_interp=filtfilt(sum(fspecial('gauss',25,10)),1,traj_interp);
% plot(tInterp,traj_interp_raw,'xr',tInterp,traj_interp)

%% create the scoring function with log scale normalized parameter inputs
input=BennettModelInput(1); 
input(17:19)=[]; 
input=input(:)'; 
fScore=@(norminput) Bennett_Optimize_Ca(10.^(norminput).*input,L,tInterp,traj_interp);

%% defint initial randomization and selection functions
Eps=1.25;
Iter=50; 
PopSize=24; 

fRand = @(N) (rand(N,32)-0.5)*4; 
% Prior = makedist('Normal','mu',0,'sigma',1); 
% fRand = @(N) random(Prior,N,32); 
fEval=@(Score,~) Score<Eps & max(abs(Pop))<4; 

%% Find initial population
[FirstGenPop,FirstGenScores,FirstGenSim] = rejectionSampling(fScore,fRand,fEval,PopSize);

%% plot first gen
figure(1)
clf
for i=1:PopSize
    subplot(5,5,i)
    try
        plot(tInterp,traj_interp,'xb',tInterp,FirstGenSim{i},'r')
        title(sprintf('%0.2f',FirstGenScores(i)))
    catch
    end
end

%% repeat for Iter iterations
AllPops=cell(Iter+1,1); 
AllPops{1}=FirstGenPop; 
AllScores=cell(size(AllPops)); 
AllScores{1}=FirstGenScores; 
Weights=cell(size(AllPops)); 
Weights{1}=ones(1,PopSize)/PopSize; 
AllSim=cell(size(AllPops)); 
AllSim{1}=FirstGenSim;


%%
t0=now; 
for i=1:Iter

    %% update the randomization and selection criteria
    
    % new random number is based on previous generation with weights
    
    % with adaptive kernbel
    fRand = @(N) AllPops{i}(randsample(PopSize,N,true,Weights{i}),:)+mvnrnd(zeros(1,32),cov(AllPops{i}),N); 
%     fRand = @(N)
%     AllPops{i}(randsample(PopSize,PopSize,true,Weights{i}),:)+randn(N,32)/2; % non adatpibe
%     fRand = @(N) AllPops{i}(randi(PopSize,N,1),:)+randn(N,32)/2;  % without weihgs
    
    % selection uses a smaller Eps 
    fEval=@(Score,Pop) Score<prctile(AllScores{i},90) & max(abs(Pop),[],2)<4;  
    
    %% call the new sampling routine
    [AllPops{i+1},AllScores{i+1},AllSim{i+1}] = rejectionSampling(fScore,fRand,fEval,PopSize);
    
    %% calculate weights 
    Weights{i+1}=nan(1,PopSize); 
    for j=1:PopSize
        Weights{i+1}(j)=Weights{i}*((mvnpdf(AllPops{i}-repmat(AllPops{i+1}(j,:),PopSize,1))).^(1/32));
        %         prod(Prior.pdf(AllPops{i+1}(j,:)))*
    end
    Weights{i+1}(max(abs(AllPops{i}),[],2)>4)=0; 
    if max(Weights{i+1})==0, 
        Weights{i+1}(:)=1/PopSize; 
    end
    %     Weights{i+1}=Weights{i+1}+eps;
    Weights{i+1}=Weights{i+1}/sum(Weights{i+1}); 
    
    %% update 
    fprintf('Finished an outter itteration: %g time: %s - yay! max: %0.2f min: %0.2f non-zeros weights: %g\n',...
            i,datestr(now-t0,13),max(mean(AllPops{i})),min(mean(AllPops{i})),nnz(Weights{i+1}));
    
end

%% plot last gen
figure(2)
for i=1:PopSize
    subplot(5,5,i)
    try
        [f,sm,f0,f1]=fScore(AllPops{end}(i,:)); 
        f0=f0/20000; 
        f1=f1/4; 
        plot(tInterp(1:10:end),traj_interp(1:10:end),'xb',tInterp,sm,'r')
        title(sprintf('%0.2f,%0.2f,%0.2f',f,f0,f1))
    catch
    end
    i
end

%% show score convergence
figure(3)
gboxplot(AllScores)

%% show weights
figure(4)
W=cat(1,Weights{:});
imagesc(W)

%% few interactive diagnostics: 
% corr of cells
figure(5),
for i=1:Iter, 
    imagesc(corr(AllPops{i}'),[-1 1]),
    title(num2str(i)),
    pause(0.33),
end
% varaibility of identified parameters
figure(6), 
for i=1:Iter, 
    boxplot(AllPops{i}),
    title(num2str(i)),
    pause(0.33),
end

%% parameters over iterations
M=cellfun(@mean,AllPops,'uniformoutput',0);
M=cat(1,M{:});
S=cellfun(@std,AllPops,'uniformoutput',0);
S=cat(1,S{:});
figure(7)
subplot(1,2,1)
imagesc(M,[-2 2])

subplot(1,2,2)
imagesc(S,[0 2])

%%
AllModels=cat(1,AllPops{:}); 
GoodModels = AllModels(cat(1,AllScores{:})<0.1,:); 

