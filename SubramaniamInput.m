function [K,x0,L,t0,t,dt] = SubramaniamInput( )
    % the vector of parameters
    K = zeros(79,1);
    % vector of initial conditions
    x0 = zeros(19,1);
    % ligand(L)
    x0(1) = 0;
    % receptor(R)
    x0(2) = 4.11e-2;
    % ligand-receptor complex(L.R)
    x0(3) = 0;
    % G protein beta gamma
    x0(4) = 8.28e-3;
    % G protein coupled receptor kinase(GRK)
    x0(5) = 2.59e-3;
    % ligand-receptor phosphorylated(L.Rp)
    x0(6) = 0;
    % receptor-phosphorylated (Rp)
    x0(7) = 0;
    % ligand-receptor internalized (L.Ri)
    x0(8) = 0;
    % receptor phosphorylated internalized (Rpi)
    x0(9) = 0;
    % receptor pool (Rpool)
    x0(10) = 0.001;
    % G protein alpha internalized - GTP complex (GaiT)
    x0(11) = 0;
    % G protein alpha internalized - GDP complex (GaiD)
    x0(12) = 8.12e-3;
    % PIP2  
    x0(13) = 110;
    % IP3
    x0(14) = 0.1;
    % CaM
    x0(15) = 3.98;
    % cytosolic calcium (Cai)
    x0(16) = .05;
    % ER calcium (Caer)
    x0(17) = 450;
    % fraction of un-inhibited ip3 receptor (h)
    x0(18) = 0.8;
    % mitochondria calcium (Camit)
    x0(19) = 0;
    %Prtot_e
    K(1) = 1.91e2;
    %Km_e
    K(2) = 2.43;
    %Prtot_x
    K(3) = 8.77;
    %Km_x
    K(4) = 1.39e-1;
    %Prtot_er
    K(5) = 6.02e4;
    %Km_er
    K(6) = 5.46e2;
    %rho_er
    K(7) = 5.33e-2;
    %kon 
    K(8) = 1.04e-1;
    %Kinh
    K(9) = 1;
    %vmax_ch 
    K(10) = 1.89e-1;
    %Kip3 
    K(11) = 1.36e-1;
    %d3 
    K(12) = 1.05;
    %Kact
    K(13) = 8.14e-2;
    %Vmax 
    K(14) = 1.14e2;
    %Kp 
    K(15) = 7.54e-1;
    %ker_leak
    K(16) = 2.03e-3;
    %Vmax_pmca_l
    K(17) = 8.93e-2;
    %Km_pmca_l
    K(18) = 1.13e-1;
    %Vmax_pmca_h
    K(19) = 5.9e-1;
    %Km_pmca_h 
    K(20) = 4.42e-1;
    %Vmax_ncx 
    K(21) = 1e-1;
    %Km_ncx 
    K(22) = 1;
    %kin 
    K(23) = 5.06e2;
    %kout
    K(24) = 4.76e2;
    %K2 
    K(25) = 9.92e-1;
    %km 
    K(26) = 1.5e-3;
    %PLCBtot 
    K(27) = 5.61e-3;
    %PIP2tot 
    K(28) = 1.29e2;
    %Xpip2_gen
    K(29) = 1.01e-1;
    %GRKtot = x0(5);
    K(30) = x0(5);
    %Km_grk = 5.07e-3;
    K(31) = 5.07e-3;
    %CaMtot = x0(15);
    K(32) = x0(15);
    %Km_cai2_cam
    K(33) = 8.36e-1;
    %Vmax_pm_ip3_dep
    K(34) = 2.26e-1;
    %kf1 
    K(35) = 5.89e1;
    %kb1
    K(36) = 2.02e-1;
    %kb2
    K(37) = 2.52e1;
    %kf3 
    K(38) = 1.22e-1;
    %Km_cai_3 
    K(39) = 1.01e-1;
    %kf4
    K(40) = 1.26e2;
    %kf5
    K(41) = 2.05e-3;
    %kf6 
    K(42) = 6.62e-4;
    %kf7
    K(43) = 1.67;
    %kf8 
    K(44) = 5.00e-3;
    %kf9
    K(45) = 0.001;
    %kf10 
    K(46) = 0.00001;
    %kf11 
    K(47) = 0.001;
    %kf12 
    K(48) = 3.67e-7;
    %kf13
    K(49) = 5e-2;
    %kf14 
    K(50) = 1.59e4;
    %kb14 
    K(51) = 4.18e-2;
    %kf15 
    K(52) = 1.38e-3;
    %Km_gid_15 
    K(53) = 6.81e-2;
    %kcat_16 
    K(54) = 1.54;
    %Km_16
    K(55) = 2.35e-1;
    %kf17
    K(56) = 5e-4;
    %kf18
    K(57) = 5.07e2;
    %Km_cai_18
    K(58) = 3.64e-1;
    %Km_gby_18 
    K(59) = 3.64e-1;
    %kf19 
    K(60) = 1.49;
    %kf20 
    K(61) = 1;
    %kf22 
    K(62) = 1.1e-2;
    %kb22 
    K(63) = 2.75e-3;
    %A 
    K(64) = 1.02e-2;
    %Arr 
    K(65) = 8.63e-3;
    %rho_m 
    K(66) = 0.01;
    %Gby_tot
    K(67) = x0(4);
    %K3 = 5;
    K(68) = 5;
    %Beta_m 
    K(69) = 0.0025;
    %kf2 = kb2/Km_grk;
    K(70) = K(37)/K(31);
    %kb5 = kf1;
    K(71) = K(35); 
    %kb21 
    K(72) = 72;
    %kf21 = kb21/Km_cai2_cam is another implementation;
    %K(73) = 86;
    K(73) = K(72)/K(33);
    %Km_pm_ip3_dep 
    K(74) = 1;
    %Vmax_pm_ip3_dep 
    K(75) = 0.2;
    %vpm_leak 
    K(76) = 0.03;
    % state variables that are calculated algebraically 
    %T 
    K(77) = 468;
    %Gaitot 
    K(78) = x0(11) + x0(12);
    %Gbytot 
    K(79) = x0(4);
    %x0(19) = x0(16)^4/((K2^4 + x0(16)^4)*(km + (kout*x0(16)^2)/(K3^2 + x0(16)^2)));
    
    
    %L is the ligand amount
    L = 0.03;
    % t0 is the period before perturbation
    t0 = 1000;
    % t is the period after perturbation
    t =  2000;
    
    % dt is the interplocated timestep
    dt = 1; 
    
end

