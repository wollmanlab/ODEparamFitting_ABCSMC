%% this function calculates the equilibrium level before stimulation
function [x0] = BennettModelEquilibrium(K)

    Cabas = K(31)/(K(34)*K(27)/(K(33)*K(28)*K(35)) -1)^(1/2);
    % Gbas, basal number of G protein molecules
    Gbas = K(10)*K(14)*K(8)/(K(10)*K(14) + K(11));
    %rhbas is the basal level of rh
    rhbas = K(16)*(Cabas/(K(15) + Cabas))*Gbas;
    % PIP2bas, basal number of PIP2 molecules
    PIP2bas = K(5)*K(13)*K(12)/(K(5)*(rhbas + K(13)) + rhbas*K(13));
    % IP3bas, basal IP3 concentration
    IP3bas = rhbas*K(13)*K(12)/(K(17)*K(18)*(K(5)*(rhbas + K(13)) + rhbas*K(13)));
    % hbas
    hbas =K(21)*(K(20) + IP3bas)/(K(20)*K(21) + IP3bas*(Cabas + K(21)) + K(22)*Cabas);
        
    
    x0 = zeros(7,1);
    % x0(1) = initial concentration for unphosphorylated  surface receptors
    %x0(1) = K(7)*K(1);
    x0(1) = K(1);
    % x0(2) = initial concentration for phosphorylated surface receptors
    x0(2) = 0;
    % x0(3) = initial concentration of activated G protein
    x0(3) = Gbas;
    %x0(4) = PIP2 initial concentration
    x0(4) =PIP2bas;
    % x0(5) = IP3 initial concentration
    x0(5) =IP3bas;
    % x0(6) = cytosolic calcuim initial concentration
    x0(6) =Cabas;
    % initial concentration of h
    x0(7) = hbas;    
end