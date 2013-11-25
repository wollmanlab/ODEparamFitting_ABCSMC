    %% MODEL TESTING => CONSERVATION
    % This script checks the validity of the Subramaniam Model by checking the conservations of the species.
    % It checks for the conservation principlles

function SubramaniamConservation(K,x0,L,t0,t,dt)

    %% Create model inputs (using the Subramniam Input function) 
    % Inputs:
    %     K = the vector of parameters
    %     x0 = initial conditions before ligand perturbation
    %     L = ligand amount for perturbation
    %     t0 = the time it takes for system to settle to steady state with zero input
    %     t = the time the system will go after t0
    % Outputs:
    %     T = vector of time points
    %     X = matrix of species concentrations

    %% Test #1 - Ligand conservation

    % Turn off the degradation flux of ligand
    KL = K;
    KL(44) = 0;
    % first get the zero ligand input condition
    [T0L,X0L,] = SubramaniamModel(KL,t0,x0,dt);
    xL = X0L(end,:);
    xL(1) = L;
    [TL,XL] = SubramaniamModel(KL,t,xL,dt);

    % plot 
    figure(1);
    clf
    subplot(1,2,1)
    plot(T0L,X0L(:,1) + X0L(:,3) + X0L(:,6) + X0L(:,8),'b',TL+ max(T0L),XL(:,1) + XL(:,3) + XL(:,6) + XL(:,8),'r');
    %axis([0 3000 -0.005 0.015])
    title('Simulation results')
    xlabel('Time (seconds)');
    ylabel('Concentration (\muM)');

    subplot(1,2,2)
    plot(T0L,zeros(size(T0L)),'b--',TL+ max(T0L),L*ones(size(TL)),'r--'); 
    %axis([0 3000 -0.005 0.015])
    title('Expected')
    suptitle('Tot. Ligand');
    xlabel('Time (seconds)');
    ylabel('Concentration (\muM)');
    h=legend('Before','After');
    set(h,'location','best'); 

    %% Test #2 - Receptor conservation
    % Check for conservation of receptor


    % Turn off the degradation flux of receptor
    KR = K;
    KR(46) = 0;

    [T0R,X0R] = SubramaniamModel(KR,t0,x0,dt);
    xR = X0R(end,:);
    xR(1) = L;
    [TR,XR] = SubramaniamModel(KR,t,xR,dt);

    figure(2)
    clf
    subplot(1,2,1);
    plot(T0R, (X0R(:,2)+ X0R(:,3) + X0R(:,6) + X0R(:,7) + X0R(:,8) + X0R(:,9) + X0R(:,10))./(X0R(1,2)+ X0R(1,3) + X0R(1,6) + X0R(1,7) + X0R(1,8) + X0R(1,9) + X0R(1,10)),'b',TR+ max(T0R),(XR(:,2)+ XR(:,3) + XR(:,6) + XR(:,7) + XR(:,8) + XR(:,9) + XR(:,10))./(X0R(1,2)+ X0R(1,3) + X0R(1,6) + X0R(1,7) + X0R(1,8) + X0R(1,9) + X0R(1,10)),'r');
    xlabel('Time (seconds)');
    ylabel('Fraction changed');
    axis([0 3000 0.999 1.001])
    title('simulation')

    subplot(1,2,2)
    plot(T0R,ones(size(T0R)),'b--',TR,ones(size(TR)),'r--')
    axis([0 3000 0.999 1.001])
    xlabel('Time (seconds)');
    ylabel('Fraction changed');
    title('expected')

    suptitle('Tot. Receptor');

    %% Test conservation of Ga 

 
    % first let the species go to steady states
    % This steady state output can be used for all the subsequent 
    [T0,X0,speciesArray,fluxMatrix0] = SubramaniamModel(K,t0,x0,dt);
    % Check for conservation of G protein alpha
    x1 = X0(end,:);
    x1(1) = L;
    [T1,X1,speciesArray1,fluxMatrix1] = SubramaniamModel(K,t,x1,dt);
    % use interpolation to get time series values in speciesArray that
    % correspond to the values in the output X vector
    GiDInterp0 = speciesArray(:,2);
    GiDInterp = speciesArray1(:,2);
    % plot the conservation for G protien alpha

    % speciesArray = [IP3p,GiD,GRK_Gby, Cai2CaMGRK, Cai2CaM];

    figure(3)
    subplot(2,4,2);
    plot(T0, (X0(:,11) + X0(:,12) + GiDInterp0)./(X0(1,11) + X0(1,12) + GiDInterp0(1)) ,'b',T1 + max(T0),(X1(:,11) + X1(:,12) + GiDInterp)./(X0(1,11) + X0(1,12) + GiDInterp0(1)),'r');
    title('Tot. G alpha');
    xlabel('Time (seconds)');
    ylabel('Fraction changed');


    % Check for conservation of G protein beta gamma
    GRKGbyInterp0 = speciesArray(:,3);
    GRKGbyInterp = speciesArray1(:,3);
    subplot(2,4,4);
    plot(T0,(X0(:,4) + GiDInterp0 + GRKGbyInterp0)./(X0(1,4) + GiDInterp0(1) + GRKGbyInterp0(1)),'b', T1 + max(T0), (X1(:,4) + GiDInterp + GRKGbyInterp)./(X0(1,4) + GiDInterp0(1) + GRKGbyInterp0(1)),'r');
    title('Tot. G beta gamma');
    xlabel('Time (seconds)');
    ylabel('Fraction changed)');

    % Check for conservation of GRK
    Cai2CaMGRKInterp0 = speciesArray(:,4);
    Cai2CaMGRKInterp  = speciesArray1(:,4);
    subplot(2,4,5);
    plot(T0,(X0(:,5) + GRKGbyInterp0 + Cai2CaMGRKInterp0)./(X0(1,5) + GRKGbyInterp0(1) + Cai2CaMGRKInterp0(1)),'b', T1 + max(T0), (X1(:,5) + GRKGbyInterp + Cai2CaMGRKInterp)./(X0(1,5) + GRKGbyInterp0(1) + Cai2CaMGRKInterp0(1)), 'r');
    title('Tot. GRK');
    xlabel('Time (seconds)');
    ylabel('Fraction changed');
    % Check for conservation of CaM
    Cai2CaMInterp0 = speciesArray(:,5); 
    Cai2CaMInterp = speciesArray1(:,5);

    subplot(2,4,6);
    plot(T0, (X0(:,15) + Cai2CaMGRKInterp0 + Cai2CaMInterp0)./(X0(1,15) + Cai2CaMGRKInterp0(1) + Cai2CaMInterp0(1)),'b',T1 + max(T0), (X1(:,15) + Cai2CaMGRKInterp + Cai2CaMInterp)./(X0(1,15) + Cai2CaMGRKInterp0(1) + Cai2CaMInterp0(1)), 'r');
    title('Tot. CaM');
    xlabel('Time (seconds)');
    ylabel('Fraction changed');
    % Check for conservation of PIP2
    ip3pInterp0 = speciesArray(:,1); 
    ip3pInterp = speciesArray1(:,1);

    subplot(2,4,7);
    plot(T0,(X0(:,13) + ip3pInterp0 + X0(:,14))./(X0(1,13) + ip3pInterp0(1) + X0(1,14)), 'b',T1 + max(T0),(X1(:,13) + ip3pInterp + X1(:,14))./(X0(1,13) + ip3pInterp0(1) + X0(1,14)), 'r');
    title('Tot. PIP2');
    xlabel('Time (seconds)');
    ylabel('Fraction changed');

    % The second plot panel has three plots 
    figure(5);
    % The second plot panel has plots on time profiles of species 
    % Plot on ligand
    subplot(3,3,1);
    plot(T0,X0(:,1),'b',T1 + max(T0),X1(:,1),'r');
    title('Ligand vs time');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on cytosolic calcium
    subplot(3,3,2);
    plot(T0,X0(:,16),'b',T1 + max(T0),X1(:,16),'r');
    title('Cytosolic calcium vs time');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on ER calcium 
    subplot(3,3,3);
    plot(T0,X0(:,17),'b',T1 + max(T0),X1(:,17),'r');
    title('ER calcium vs time');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on mitochondria calcium
    subplot(3,3,4);
    plot(T0,X0(:,19),'b',T1 + max(T0),X1(:,19),'r');
    title('mitochondria calcium vs time');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on IP3 
    subplot(3,3,5);
    plot(T0,X0(:,19),'b',T1 + max(T0),X1(:,19),'r');
    title('mitochondria calcium vs time');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on G protein beta gamma
    subplot(3,3,6);
    plot(T0,X0(:,4),'b',T1 + max(T0),X1(:,4),'r');
    title('G protein beta gamma vs time');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % The third plot panel has plots on fluxes
    figure(6);
    % Plot on IP3 generation flux
    subplot(2,2,1);
    plot(T0,fluxMatrix0(1,:),'b',T1 + max(T0),fluxMatrix1(1,:),'r');
    title('IP3 Flux');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jch, the ip3 channel calcium flux
    subplot(2,2,2);
    plot(T0,fluxMatrix0(2,:),'b',T1 + max(T0),fluxMatrix1(2,:),'r');
    title('Jch');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jserca
    subplot(2,2,3);
    plot(T0,fluxMatrix0(3,:),'b',T1 + max(T0),fluxMatrix1(3,:),'r');
    title('Jserca');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jer,leak
    subplot(2,2,4);
    plot(T0,fluxMatrix0(4,:),'b',T1 + max(T0),fluxMatrix1(4,:),'r');
    title('Jer-leak');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jpm,ip3dep
    figure(7);
    subplot(2,2,1);
    plot(T0,fluxMatrix0(5,:),'b',T1 + max(T0),fluxMatrix1(5,:),'r');
    title('Jpm-ip3dep');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jpm,leak
    subplot(2,2,2);
    plot(T0,fluxMatrix0(6,:),'b',T1 + max(T0),fluxMatrix1(6,:),'r');
    title('Jpm-leak');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jpmca
    subplot(2,2,3);
    plot(T0,fluxMatrix0(7,:),'b',T1 + max(T0),fluxMatrix1(7,:),'r');
    title('Jpmca');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jncx
    subplot(2,2,4);
    plot(T0,fluxMatrix0(8,:),'b',T1 + max(T0),fluxMatrix1(8,:),'r');
    title('Jncx');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jmit,in 
    figure(8);
    subplot(1,2,1);
    plot(T0,fluxMatrix0(9,:),'b',T1 + max(T0),fluxMatrix1(9,:),'r');
    title('Jmit-in');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');
    % Plot on Jmit,out
    subplot(1,2,2);
    plot(T0,fluxMatrix0(10,:),'b',T1 + max(T0),fluxMatrix1(10,:),'r');
    title('Jmit-out');
    xlabel('Time (seconds)');
    ylabel('Concentration(\muM)');

 
end

