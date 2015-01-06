%% This script performs the trajectory sensitivity analysis by approximating the hessian matrix
%% Init
clear 
close all
clc
%% load experimental data
S= load ('/home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/traj_147.mat')
traj_interp = S.traj_interp;
S= load('/home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/tInterp.mat')
tInterp = S.tInterp;
L=10; 

%% Load the optimal fits for trajectory # 147
optimizeFolderPath = '/home/jason/calciummodeling/Traj_Compare_Conditions/';
folderDir = '/home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/';
resultsOrigCaArray = cell(1,1);
count=1;
origCaIndArray = zeros(1,1);
for i =1:595
    if exist(strcat(optimizeFolderPath,'Traj_Orig_Ca/results_orig_ca_',num2str(i),'.mat'),'file')
        S = load(strcat(optimizeFolderPath,'Traj_Orig_Ca/results_orig_ca_',num2str(i),'.mat'));
        resultsOrigCaArray{count} = S.Results;
        origCaIndArray(count) = i;
        count = count +1;
    end
end
paramMat = zeros(length(resultsOrigCaArray),length(resultsOrigCaArray{1}.xbest));
for i = 1:size(paramMat,1)
    resultsVec = resultsOrigCaArray{i};
    paramMat(i,:) = resultsVec.xbest;
end
%% initialize the cell array for Hessian matrix
hessianArray = cell(size(paramMat,1),1);

%% decide the parameter set to use for sensitivity analysis
differential = 0.1;
t = max(tInterp);
dt = tInterp(2) - tInterp(1);
for i =1:length(hessianArray)
    disp(strcat('i = ',num2str(i)));
    input = paramMat(i,:);
    ver =1;
    [K,~,~,~,~] = BennettModelInput(ver);
    K(1:16) = input(1:16);
    K(20:end) = input(17:end);
    x0 = BennettModelEquilibrium(K);
    if ~isreal(x0)
        [~,~,~,~,x00] = BennettModelInput(ver);
        t0 = 5000;
        [T0,X0] =BennettModelSBML(K,t0,dt,0,x00);
        x0 = X0(end,:);
    else
        t0 = 1000;
        [T0,X0] =BennettModelSBML(K,t0,dt,0,x0);
        x0 = X0(end,:);
    end
    [jac] = BennettJacobian(differential,K,t0,t,dt,L,x0); 
    hessian = jac'*jac;
    hessianArray{i} = hessian;
    
end

%% Decompose the Hessian matrix into its eigenvectors and eigenvalues
% the cell array is an array of structs. Each of the structs contains the
% sorted eigenvalues, three most significant paramaters in the respective
% eigenvectors, and the norm of their values in eigenvectors. 
matCell = cell(length(hessianArray),1);
for i =1:length(hessianArray)
    hessian = hessianArray{i};
    [V,D] = eig(hessian);
    % flip the V and D left-right
    V = fliplr(V);
    D = fliplr(D);
    % declare fileds of struct
    field1 = 'eigenvalues';
    field2 = 'mostSigParam';
    field3 = 'norm';
    values1 = sum(D,1);
    values2 = zeros(3,size(V,2));
    values3 = zeros(3,size(V,2));
    for j =1:size(values2,2)
        eigenVec = V(j,:);
        absEigenVec = abs(eigenVec);
        [sortedEigenVec,I] = sort(absEigenVec,'descend'); 
        % Get the  
        values2(:,j) = I(1:3);
        values3(:,j) = sortedEigenVec(1:3);
    end
    hessianStruct = struct(field1,values1,field2,values2,field3,values3);
    matCell{i} = hessianStruct;
end

