% This optimization problem only varies the number of proteins plus the
% kinetic parameters in the receptor regulation module, using the original
% Bennett model
function [F,caSim,F0,F1] = Bennett_Optimize_Ca(input,L,tInterp,traj_interp)



    % read in the input parameters to prepare for simulation
    ver =1;
    [K,~,~,~,x0] = BennettModelInput(ver);
    K(1:16) = input(1:16);
    K(20:end) = input(17:end);
    
    K(7)=min(K(7),1);
    K(7)=max(K(7),0);
%     x0 = BennettModelEquilibrium(K);
    t0 = 5000;
    t = max(tInterp);
    dt = tInterp(2) -tInterp(1);
    % carry the simulation to equilibrium
    %tEqui = [0:dt:t0]';
    [T0,X0,runwasok] = BennettModel(K,t0,dt,0,x0);
    if ~runwasok
        F=inf; 
        caSim=nan; 
        F0=inf; 
        F1=inf; 
        return
    end
    x01 = X0(end,:)';
    % execute the simulation
    %tVecSim = [0:dt:t]';
    [T,X,runwasok] = BennettModel(K,t,dt,L,x01);
    if ~runwasok
        F=inf;
        caSim=nan;
        F0=inf;
        F1=inf;
        return
    end
    % assemble the time vector and the calcium trace from the simulation
    caSim = X(:,6)./(X(1,6));
    % get the differentials of the vectors and compare between the two
    diffTraj = diff(traj_interp);
    diffSim = diff(caSim); 
    % compare the least squares of the two vectors and also the least
    % squares of time points
%     if length(diffTraj) ~= length(diffSim)
%         F0 = inf;
%         F1 = inf;
%         
%     else
        diffVec = (diffSim - diffTraj).^2;
        %diffVec = ((diffSim - diffTraj)./diffSim).^2;
        % weight is the multiplier to the first 1/3 of the vector
        weight = 1;
        diffVec(1:floor(length(diffVec)/3)) = diffVec(1:floor(length(diffVec)/3))*weight;
        F1 = sum(diffVec)/6; 
        F0 = sum((traj_interp-caSim).^2)/20000;
%     end
    
    F=F0+F1;
% F=F0; 
   
    
%     % Equality constraints. declared before inequality constraints
%     % no negative entries in the time series
%     g(1) = length(find(X(:,1) <0) );
%     g(2) = length(find(X(:,2) <0) );
%     g(3) = length(find(X(:,3) <0)  );
%     g(4) = length(find(X(:,4) <0) );
%     g(5) = length(find(X(:,5) <0) );
%     g(6) = length(find(X(:,6) <0) );
%     g(7) = length(find(X(:,7) <0)  );
%     % no negative entries in the time sereis as simlation goes to
%     % equilibrium
%     g(8) = length(find(X0(:,1) <0) );
%     g(9) = length(find(X0(:,2) <0));
%     g(10) = length(find(X0(:,3) <0) );[
%     g(11) = length(find(X0(:,4) <0));
%     g(12) = length(find(X0(:,5) <0) );
%     g(13) = length(find(X0(:,6) <0));
%     g(14) = length(find(X0(:,7) <0));
%     g(15) = x01(3);
%     g(16) = x01(6);
%     g(17) = x01(7);
    
    
    
return
