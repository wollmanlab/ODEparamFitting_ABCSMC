%% Init
clear
close all
clc

% %% Calculate few critical expressions
% syms Cmin Mw a t r h D S positive
% rmn = solve(Cmin == Mw*a/h/4/pi/D/t*exp(-r^2/(4*D*t)),r);
% r1=rmn(1); 
% pretty(simplify(r1))
% dr1dt = diff(r1,t);
% t1mx = solve(dr1dt,t);
% pretty(simplify(t1mx));
% r1mx = subs(r1,t,t1mx);
% pretty(simplify(r1mx))
% 
% Mr1 = simplify(S*pi*r1mx^2);
% pretty(simplify(Mr1))
% 
% %% with density dependent degradation
% syms Kdeg positive
% rmn = solve(Cmin == exp(-S*Kdeg*t) * Mw*a/h/4/pi/D/t*exp(-r^2/(4*D*t)), r); 
% r2 = rmn(1); 
% matlabFunction(r2)

%% do some substitutaitons

% 
% % find depndnecy over Mw - const 
% figure(1)
% clf
% MwVec = 20:10:350;
% 
% KdegVec=[0 5 100];
% MrVec=zeros(numel(MwVec),numel(KdegVec)); 
% for j=1:numel(KdegVec)
%     Kdeg=KdegVec(j);
%     for i=1:numel(MwVec)
%         Mw = MwVec(i);
%         r = real(sqrt(4*D*Tvec.*(log(Mw*a./h./Cmin/4/pi/D./Tvec)-S*Kdeg.*Tvec)));
%         r(r==0)=nan;
%         MrVec(i,j)=pi*S*max(r)^2;
%     end
% end
% plot(MwVec,MrVec)
% xlabel('# Wounded cells (Density is constant => changing wound size)')
% ylabel('# responding cells')
% h=legend('No deg (Kdeg=0)','Some deg (Kdeg=5)','No deg (KDeg = 100'); 
% set(h,'location','best')

%%
D=200; % ATP diffusion
Cmin = 600; % molecules per micron cube that will result in 10% of cell response. 
% Cmin = 600 => 0.1 uM
h = 100; % height of cylinder

% a = ~1E9 molecules per cell
% Similar calculation could be done based on internal concentration of ATP
% within a cell 
% based on (bionumbers) ~1 mM concentration (1E6 nM)
%                       ~1000 molecules per cell per 1nM



% S = 1/750; % cells / um^2

Tvec = 0:180000;
figure(2)
set(2,'position',[680   164   890   934])
clf
a=5E8; % molecules per cell (~1mM)
Kdeg=0;
TotArea = (6.5/7)^2*2048*2065*0.8; 
SVec = 1./linspace(500,2000,100); 
MwArea = pi*([150 300 450]/2).^2;
MrVec=nan(numel(SVec),numel(MwArea));

for j=1:numel(SVec)
    S=SVec(j);
    for i=1:numel(MwArea)
        Mw = MwArea(i)*S;
        r = real(sqrt(4*D*Tvec.*(log(Mw*a./h./Cmin/4/pi/D./Tvec)-S*Kdeg.*Tvec)));
        r(r==0)=nan;
        [rmx,ix]=max(r);
        Tmx0(j,i)=Tvec(ix); 
        MrVec(j,i)=pi*S*rmx^2;
    end
end
subplot(3,2,1)
MwMat = repmat(MwArea,numel(SVec),1).*repmat(SVec(:),1,numel(MwArea)); 
plot(MwMat,MrVec)
xlabel('# cell crushed (Changing Density, i.e. response^2)')
ylabel('# responding cells')
title('No Degradation')
h2=legend('150','300','450'); 
set(h2,'location','northwest')


subplot(3,2,3)
PerResponse = MrVec./(TotArea*repmat(SVec(:),1,numel(MwArea)));
plot(MwMat,PerResponse)
xlabel('# cell crushed (Changing Density, i.e. response^2)')
ylabel('% of cell responding responding cells')
title('No Degradation')
ylim([0 1])


Kdeg=2;
fold=4; 
a=a*fold; 
MrVec=zeros(numel(SVec),numel(MwArea)); 
Tmx = zeros(size(MrVec)); 
for j=1:numel(SVec)
    S=SVec(j);
    for i=1:numel(MwArea)
        Mw = MwArea(i)*S;
        r = real(sqrt(4*D*Tvec.*(log(Mw*a./h./Cmin/4/pi/D./Tvec)-S*Kdeg.*Tvec)));
        r(r==0)=nan;
        [rmx,ix]=max(r);
        TmxDeg(j,i)=Tvec(ix); 
        MrVec(j,i)=pi*S*rmx^2;
    end
end
subplot(3,2,2)
plot(MwMat,MrVec)
xlabel('# cell crushed (Changing Density, i.e. response^2)')
ylabel('# responding cells')
title(sprintf('Kdeg = %g [mol/cell/sec] & Conc = * %g-fold',Kdeg,fold))

subplot(3,2,4)
PerResponse = MrVec./(TotArea*repmat(SVec(:),1,numel(MwArea)));
plot(MwMat,PerResponse)
xlabel('# cell crushed (Changing Density, i.e. response^2)')
ylabel('% of cell responding responding cells')
title(sprintf('Kdeg = %g [mol/cell/sec] & Conc = * %g-fold',Kdeg,fold))
ylim([0 1])

% look at timing
subplot(3,2,5)
plot(MwMat,Tmx0)
axis([0 350 0 600])
xlabel('# cell crushed (Changing Density, i.e. response^2)')
ylabel('Time till the last responds')
title('No deg')

subplot(3,2,6)
plot(MwMat,TmxDeg)
title(sprintf('Kdeg = %g [mol/cell/sec] & Conc = * %g-fold',Kdeg,fold))
axis([0 350 0 600])
xlabel('# cell crushed (Changing Density, i.e. response^2)')
ylabel('Time till the last responds')

%% start substituting values... 

% dr2dt = diff(r2,t); 
% t2mx = solve(dr2dt,t); 
% disp('t2mx: ');
% pretty(simplify(t2mx))
% r2mx = subs(r2,t,t2mx); 
% disp('r2mx: ');
% pretty(simplify(r2mx))
% 
% Mr2 = S*pi*r2mx^2; 
% pretty(simplify(Mr2))