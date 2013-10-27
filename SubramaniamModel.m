%{
This function simulates the Subramaniam et al Ca2+ model

Inputs:
    K = the vector of parameters
    x0 = vector of initial conditions
    x0(1) = ligand = varied form 1nM to 500 nM
    x0(2) = Receptor (R) = 4.11e-2
    x0(3) = ligand-receptor complex (L.R) =0
    x0(4) = G protein beta gamma = 8.28e-3
    x0(5) = G protein coupled receptor kinase (GRK) = 2.59e-3
    x0(6) = ligand-receptor-phosphorylated (L.Rp) = 0
    x0(7) = receptor-phosphorylated(Rp) = 0
    x0(8) = ligand-receptor-internalized (L.Ri) =0 
    x0(9) = receptor-phosphorylated-internalized( Rpi) = 0
    x0(10)= receptor pool(Rpool) = 0.001
    x0(11) = G protein alpha internalized -GTP complex(GaiT) = 0
    x0(12) = G protein alpha internalized -GDP complex(GaiD) = 8.12e-3
    x0(13) = PIP2 = 110
    x0(14) = IP3 = 0.1
    x0(15) = CaM = 3.98 
    x0(16) = cytosolic calcium (Cai) = <2
    x0(17) = ER calcium (Caer) = < 2000
    x0(18) = fraction of un-inhibited ip3 receptor(h) = 0.8
    x0(19) = mitochondria calcium (Camit) = 0

Outputs:


%}

function [T,X ,speciesArray,fluxMatrix] = SubramaniamModel( K,t,x0,dt)
    Prtot_e = K(1);
    Km_e = K(2);
    Prtot_x = K(3);
    Km_x = K(4);
    Prtot_er = K(5);
    Km_er = K(6);
    rho_er = K(7);
    kon = K(8);
    Kinh = K(9);
    vmax_ch = K(10);
    Kip3 = K(11);
    d3 = K(12);
    Kact = K(13);
    Vmax = K(14);
    Kp = K(15);
    ker_leak = K(16);
    Vmax_pmca_l = K(17);
    Km_pmca_l = K(18);
    Vmax_pmca_h = K(19);
    Km_pmca_h = K(20);
    Vmax_ncx  = K(21);
    Km_ncx = K(22);
    kin = K(23);
    kout = K(24);
    K2 = K(25);
    km = K(26);
    PLCBtot = K(27);
    PIP2tot = K(28);
    Xpip2_gen = K(29);
    GRKtot = K(30);
    CaMtot = K(32);
    Vmax_pm_ip3_dep = K(34);
    kf1 = K(35);
    kb1 = K(36);
    kb2 = K(37);
    kf3 = K(38);
    Km_cai_3 = K(39);
    kf4 = K(40);
    kf5 = K(41);
    kf6 = K(42);
    kf7 = K(43);
    kf8 = K(44);
    kf9 = K(45);
    kf10 = K(46);
    kf11 = K(47);
    kf12 = K(48);
    kf13 = K(49);
    kf14 = K(50);
    kb14 = K(51);
    kf15 = K(52);
    Km_gid_15 = K(53);
    kcat_16 = K(54);
    Km_16 = K(55);
    kf17 = K(56);
    kf18 = K(57);
    Km_cai_18 = K(58);
    Km_gby_18 = K(59);
    kf19 = K(60);
    kf20 = K(61);
    kf22 = K(62);
    kb22 = K(63);
    A = K(64);
    Arr = K(65);
    rho_m = K(66);
    Gby_tot = K(67);
    K3 = K(68);
    Beta_m = K(69);
    kf2 = K(70);
    kb5 = K(71);
    kb21 = K(72);
    kf21 = K(73);
    Km_pm_ip3_dep = K(74);
    Vmax_pm_ip3_dep = K(75);
    vpm_leak = K(76);
    Tk = K(77);
    Gaitot = K(78);
    Gby_tot = K(79);
    
    Tgrid = (0:dt:t)'; 
    function [dx_dt,speciesArray,fluxMatrix] = F(t,x) %#ok<INUSL>
        % Calculate the related parameters
        IP3p = PIP2tot - x(14,:) - x(13,:);
        Q = Kinh.*(x(14,:) +Kip3)./(x(14,:) + d3);
        Beta_i = (1 + Prtot_e.*Km_e./(Km_e + x(16,:)).^2 + Prtot_x.*Km_x./(Km_x + x(16,:)).^2).^-1;
        Beta_er = (1 + Prtot_er.*Km_er./(Km_er + x(17,:)).^2).^-1;
        Popen = ((x(14,:)./(x(14,:) + Kip3)).*(x(16,:)./(x(16,:) + Kact)).*x(18,:)).^3;
        GiD = Gaitot - x(11,:) -x(12,:); 
        GRK_Gby = Gby_tot - x(4,:) - GiD;
        Cai2CaMGRK = GRKtot -x(5,:) -GRK_Gby;
        Cai2CaM = CaMtot - x(15,:) - Cai2CaMGRK;
        speciesArray = [IP3p,GiD,GRK_Gby, Cai2CaMGRK, Cai2CaM,Q];
        
        %initialize the vector for fluxes
        v = zeros(31,size(x,2));
        % L + R <-> LR
        % v(1) = kf1[L][R] - kb1[LR]
        v(1,:) = kf1.*x(1,:).*x(2,:) - kb1.*x(3,:);
        % G protein-beta-gamma + GRK <-> GRK.Gprotein-beta-gamma
        % v(2) = kf2[G-beta-gamma][GRK] - kb2[GRK.G-beta-gamma]
        v(2,:) = kf2.*x(4,:).*x(5,:) - kb2.*GRK_Gby;
        %[GRK;Cai]L.R -> L.Rp
        % v3 
        v(3,:) = kf3.*x(3,:).*x(5,:).*(x(16,:)./(Km_cai_3 + x(16,:)));
        %[GRK.G_beta_gamma]LR -> L.Rp
        %v4 = kf4[LR][GRK.G_beta_gamma]
        v(4,:) = kf4.*x(3,:).*GRK_Gby;

        % L.Rp <-> L + Rp
        % v5 = kf5[L.Rp] - kb5[L][Rp]
        v(5,:) = kf5.*x(6,:) - kb5.*x(1,:).*x(7,:);

        % L.Rp -> L.Ri
        % v6
        v(6,:) = kf6.*x(6,:);
   
        % [Arr]L.Rp -> L.Ri
        %v7
        v(7,:) = kf7.*x(6,:).*Arr; 
   
        
        % L.Ri -> Li + Rpi
        % v8
        v(8,:) = kf8.*x(8,:);  
 
        %Rpi -> R
        % v9
        v(9,:) = kf9.*x(9,:); 
     
        %Rpi -> Rvac
        %v10
        v(10,:) = kf10.*x(9,:);
    
        % Rpool -> R
        % v11
        v(11,:) = kf11.*x(10,:);
        % GiD + T -> G_alpha_iT + G_beta_gamma + D
        % v12
        v(12,:) = kf12.*GiD.*Tk;
        % GaiT -> GaiD + Pi
        %v13
        v(13,:) = kf13.*x(11,:);
        %GaiD + Gby <-> GiD 
        %v14
        v(14,:) = kf14.*x(12,:).*x(4,:) - kb14.*GiD;
        % [LR]GiD + T -> GaiT + Gby + D
        %v15
        v(15,:) = kf15.*(GiD./(Km_gid_15+ GiD)).*Tk.*x(3,:);
        % [A]GaiT -> GaiD + Pi
        %v16
        v(16,:) = kcat_16.*A.*x(11,:)./(Km_16 + x(11,:));
        %PIP2 -> IP3 + DAG
        %v17
        v(17,:) = kf17.*x(13,:);
        %[Cai;Gby;PLCBtot]PIP2 -> IP3 + DAG
        % v18
        v(18,:) = kf18.*x(13,:).*PLCBtot.*(x(16,:)./(Km_cai_18 + x(16))).*(x(4,:)./(Km_gby_18 + x(4,:)));
        % IP3 -> IP3p
        %v19
        v(19,:) = kf19.*x(14,:);
        % [Xpip2 gen]IP3,p ->PIP2
        %v20
        v(20,:) = kf20.*IP3p.*Xpip2_gen;
        % 2Cai + CaM <-> Cai2.*CaM
        % v21
        v(21,:) = kf21.*x(16,:).^2.*x(15,:) - kb21.*Cai2CaM; 
        %Cai2.CaM + GRK <-> Cai2.CaM.GRK
        %v22
        v(22,:) = kf22.*Cai2CaM.*x(5,:) - kb22.*Cai2CaMGRK;
        
        
        %Jch 
        v(23,:) = vmax_ch.*Popen.*(x(17,:) - x(16,:)); 
        %Jserca
        v(24,:) = Vmax.*x(16,:).^2./(x(16,:).^2 + Kp.^2);
        %Jer,leak
        v(25,:) = ker_leak.*(x(17,:) - x(16,:));
        % Jpm,ip3dep
        v(26,:) = Vmax_pm_ip3_dep.*x(14,:).^2./(Km_pm_ip3_dep.^2 + x(14,:).^2);
        %Jpm,leak
        v(27,:) = vpm_leak;
        %Jpmca
        v(28,:) = Vmax_pmca_l.*x(16,:).^2./(x(16,:).^2 + Km_pmca_l.^2) + Vmax_pmca_h.*x(16,:).^5./(x(16,:).^5 + Km_pmca_h.^5);
        %Jncx
        v(29,:) = Vmax_ncx.*x(16,:)./(x(16,:) + Km_ncx);
        %Jmit,in
        v(30,:) = kin.*x(16,:).^4./(K2.^4 + x(16,:).^4);
        %Jmit,out
        v(31,:) = (kout.*x(16,:).^2./(K3.^2 + x(16,:).^2) + km).*x(19,:);
       
        % Initialize the dx_dt vector that keeps track of the differentials 
        dx_dt = zeros(19,size(x,2));
        
        %rate of change of ligand
        vec = zeros(1,31);
        vec(1) = -1;
        vec(5) =1;
        dLdt  = vec*v;
        dx_dt(1,:) = dLdt;
        % rate of change of receptors
        vec = zeros(1,31);
        vec(1) = -1;
        vec(9) = 1;
        vec(11)= 1;
        dRdt = vec*v;
        dx_dt(2,:) = dRdt;
        %rate of change of L.R
        vec = zeros(1,31);
        vec(1) = 1;
        vec(3) = -1;
        vec(4) = -1;
        dLRdt = vec*v;
        dx_dt(3,:) = dLRdt;
        %rate of change of Gby
        vec = zeros(1,31);
        vec(2) = -1;
        vec(12) = 1;
        vec(14) = -1;
        vec(15) = 1;
        dGbydt = vec*v;
        dx_dt(4,:) = dGbydt;
        %rate of change of GRK
        vec = zeros(1,31);
        vec(2) = -1;
        vec(22) =  -1;
        dGRKdt = vec*v;
        dx_dt(5,:) = dGRKdt;
        %rate of change of L.Rp
        vec = zeros(1,31);
        vec(3) = 1; 
        vec(4) = 1;
        vec(5) = -1;
        vec(6) = -1;
        vec(7) = -1;
        dLRpdt = vec*v;
        dx_dt(6,:) = dLRpdt;
        % rate of change of Rp
        vec = zeros(1,31);
        vec(5) = 1;
        dRpdt = vec*v;
        dx_dt(7,:) = dRpdt;
        % rate of change of L.Ri
        vec = zeros(1,31);
        vec(6) = 1;
        vec(7) = 1;
        vec(8) = -1;
        dLRidt = vec*v;
        dx_dt(8,:) = dLRidt;
        % rate of change of Rpi
        vec = zeros(1,31);
        vec(8) = 1;
        vec(9) = -1;
        vec(10) = -1;
        dRpidt = vec*v;
        dx_dt(9,:) = dRpidt;
        % rate of change of Rpool
        vec = zeros(1,31);
        vec(11) = -1;
        dRpooldt = vec*v;
        dx_dt(10,:) = dRpooldt;
        
       
        %rate of change of GaiT
        vec = zeros(1,31);
        vec(12) = 1;
        vec(13) = -1;
        vec(15) = 1;
        vec(16) = -1;
        dGaiTdt = vec*v;
        dx_dt(11,:) = dGaiTdt;
        % rate of change of GaiD
        vec = zeros(1,31);
        vec(13) = 1;
        vec(14) = -1;
        vec(16) = 1;
        dGaiDdt = vec*v;
        dx_dt(12,:) = dGaiDdt;
        %}
        
        % rate of change of PIP2
        vec = zeros(1,31);
        vec(17) = -1;
        vec(18) = -1;
        vec(20) = 1;
        dPIP2dt = vec*v;
        dx_dt(13,:) = dPIP2dt;
        %rate of change of IP3
        vec = zeros(1,31);
        vec(17) = 1;
        vec(18) = 1;
        vec(19) = -1;
        dIP3dt = vec*v;
        dx_dt(14,:) = dIP3dt;
        %rate of change of CaM
        vec = zeros(1,31);
        vec(21) = -1;
        dCaMdt = vec*v;
        dx_dt(15,:) =dCaMdt;
        
        % rate of change of Cai
        vec = zeros(1,31);
        % Jch 
        vec(23) =1;
        %Jerleak
        vec(25) = 1;
        %Jpm_ip3_dep
        vec(26) = 1;
        %Jserca
        vec(24) =-1;
        %Jpmca
        vec(28) = -1;
        %Jncx 
        vec(29) = -1;
        %Jpmleak
        vec(27) = 1;
        %Jmit_out
        vec(31) = 1;
        %Jmit_in
        vec(30) = -1;
        %v21
        vec(21) =  -2;
        dCaidt = Beta_i.*(vec*v);
        dx_dt(16,:) = dCaidt;
        % rate of change of Caer
        vec = zeros(1,31); 
        vec(24) = 1;
        vec(23) =-1;
        vec(25) = -1;
        dCaerdt = (Beta_er/rho_er).*(vec*v);
        dx_dt(17,:) = dCaerdt;
        % rate of change of h 
        dhdt = kon.*(Q -(x(16,:)+ Q).*x(18,:));
        dx_dt(18,:) = dhdt;
        % rate of change of Camit
        vec = zeros(1,31);
        vec(30) =1;
        vec(31) = -1;
        dCamitdt = (Beta_m/rho_m)*(vec*v);
        dx_dt(19,:) = dCamitdt;
        
        % Assemble the fluxMatrix which has the fluexes of IP3 generation, Jch, Jserca, Jerleak, Jpmip3dep, Jpmleak,  Jpmca, Jncx, Jmitin, Jmitout 
        fluxMatrix = [v(17,:)+v(18,:); v(23,:); v(24,:) ; v(25,:); v(26,:) ; v(27,:) ; v(28,:) ; v(29,:) ; v(30,:) ; v(31,:)];

    end
    options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Refine',1000);
    [T,X] = ode15s(@F,Tgrid,x0,options);
    [~,speciesArray,fluxMatrix] = F(Tgrid,X'); 

end

