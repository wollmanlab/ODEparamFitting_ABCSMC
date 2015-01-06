% This optimization problem only varies the number of proteins plus the
% kinetic parameters in the receptor regulation module, using the original
% Bennett model
function [F,caSim,g] = BennettABCOptimize(input,L,tInterp,traj_interp)
    % read in the input parameters to prepare for simulation
    ver =1;
    [K,~,~,~,~] = BennettModelInput(ver);
    K(1:16) = input(1:16);
    K(20:end) = input(17:end);
    
    x0 = BennettModelEquilibrium(K);
    x01 = zeros(size(x0));
    dt = tInterp(2) - tInterp(1);
    t = max(tInterp);
    if ~isreal(x0)
        [~,~,~,~,x00] = BennettModelInput(ver);
        t0 = 5000;
        [T0,X0] =BennettModelSBML(K,t0,dt,0,x00);
        x01 = X0(end,:);
    else
        t0 = 1000;
        [T0,X0] =BennettModelSBML(K,t0,dt,0,x0);
        x01 = X0(end,:);
    end
 
    % execute the simulation
    [T,X] =BennettModelSBML(K,t,dt,L,x01);
    
    % assemble the time vector and the calcium trace from the simulation
    caSim = X(:,6)./(X(1,6));
    % get the differentials of the vectors and compare between the two
    diffTraj = diff(traj_interp);
    diffSim = diff(caSim); 
    if length(diffTraj) ~= length(diffSim)
        F = inf;
    else
        diffVec = (diffSim - diffTraj).^2;
        F = sum(diffVec); 
    end 
   
    
%     % Equality constraints. declared before inequality constraints
%     % no negative entries in the time series
    g(1) = length(find(X(:,1) <0) );
    g(2) = length(find(X(:,2) <0) );
    g(3) = length(find(X(:,3) <0)  );
    g(4) = length(find(X(:,4) <0) );
    g(5) = length(find(X(:,5) <0) );
    g(6) = length(find(X(:,6) <0) );
    g(7) = length(find(X(:,7) <0)  );
% no negative entries in the time sereis as simlation goes to
%     % equilibrium
    g(8) = length(find(X0(:,1) <0) );
    g(9) = length(find(X0(:,2) <0));
    g(10) = length(find(X0(:,3) <0) );
    g(11) = length(find(X0(:,4) <0));
    g(12) = length(find(X0(:,5) <0) );
    g(13) = length(find(X0(:,6) <0));
    g(14) = length(find(X0(:,7) <0));
    g(15) = x01(3);
    g(16) = x01(6);
    g(17) = x01(7);
    
    
    
return
