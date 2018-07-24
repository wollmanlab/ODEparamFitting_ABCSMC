function KL = estimateKLdivergenceBasedOnNN(P)
% estimateKLdivergenceBasedOnNN calculates KL based on who the NN is for
% two samples
%   P is a cell array where each col is a sample (can have many rows)
%   KL is the calculated KL divergence 


%% estimate KL divergence between any two cell parameter samples
D=nan(numel(P)); 
N=cellfun(@(p) size(p,2),P);
for i=1:numel(P)
    D(i,i)=0.5; 
    parfor j=i+1:numel(P)
        n=min(N([i j])); 
        prm=randperm(size(P{i},2)); 
        Pi=P{i}(:,prm(1:n)); 
        prm=randperm(size(P{j},2)); 
        Pj=P{j}(:,prm(1:n)); 
        Pij=[Pi Pj]; 
        IDX=annquery(Pij,Pij,2);
        IDX=IDX(2,:); 
        D(i,j)=mean([IDX(1:size(P{i},2))<size(P{i},2) IDX(size(P{i},2):1:end)>size(P{i},2)]);  
    end
    i
end
D=1-D'; 
KL=D.*log(D*2)+(1-D).*log((1-D).*2);
for i=1:size(D,1)
    for j=2:size(KL,2)
        KL(i,j)=KL(j,i); 
    end
end

