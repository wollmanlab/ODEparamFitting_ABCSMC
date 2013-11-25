%% Init
clear
close all
clc

load P

%% Here we approximate 3D diffusion in a disc by calculating a 3D case based
% on 3D gaussian and projecting it into a disc with 100 uM height

%% Parameters
D=200; % ATP diffusion
N=1E11; % number of molecules secreted from a ~300 um diameter wound. 
Cmin = 600; % molecules per micron cube that will result in 10% of cell response. 
h = 100; % height of cylinder

% estimates of ~100 uM around the wound, assume volume of 300*300*11 um^3
% this means that we have a volume of ~1E6 um^3. 
% using approximation of ~0.6 molecules / 1um^3 = 1nM we get that number of
% molecules of ATP realease to be ~1E6*1E5 ~1E11. 

% Similar calculation could be done based on internal concentration of ATP
% within a cell 
% based on (bionumbers) ~1 mM concentration (1E6 nM)
%                       ~1000 molecules per cell per 1nM
%                       ~100 cell death
%                       6+3+2 = 11

%% grid
grd = linspace(-4000,4000,101); % in um
dx = mean(diff(grd));  
[x,y,z]=meshgrid(grd,grd,grd); 

[X,Y]=meshgrid(grd,grd); 

r=sqrt(x.^2+y.^2+z.^2);              

%% Diffusion in a cylinder  
funcC = @(r,t) N/(4*pi*D.*t.*h)*exp(-r.^2/4/D./t); 


% equation for distance that a "Front" of Cmin
% concentration diffused in. 
f_r2 = @(t) sqrt(4*D*t.*log(N./4/D./t/pi/Cmin/h)); % 100 um height, as if we divide N to 100 equal 1um height layers and diffuse... 


funcC2d = @(r,t) N/(4*pi*D.*t./h)*exp(-r.^2/4/D./t); 

%% diffusion plot in different distances and wound size
figure(123)
funcC2d = @(r,t,N) N./(4*pi.*D.*t./h).*exp(-r.^2/4/D./t)/0.6/1000;
clf, 
r0=300; 
N0=1E9; 
clr=jet(10); 
subplot(1,2,1), hold on
v = logspace(-1,1,9); 
t=0:1200; 
for i=1:numel(v)
    plot(t,funcC2d(v(i)*r0,t,N0),'color',clr(i,:)),
end
title('Distance from wound')
axis([0 max(t) 0 1000])
subplot(1,2,2), hold on
for i=1:numel(v), 
    plot(t,funcC2d(r0,t,N0./v(i)),'color',clr(i,:)),
end
title('Wound size')
axis([0 max(t) 0 1000])

%% calculate C on the grid r
dst = sqrt(X.^2+Y.^2);
[dst,ix] = unique(dst(:));
dgrid = linspace(0,1500,351); 
Tgrid = linspace(0,1200,1201); 
% Ctime = zeros(numel(Tgrid),numel(dgrid)); 
% for i=1:numel(Tgrid)
%     C = funcC(r,Tgrid(i)); % is in number of molecules 
%     Cdisc = sum(C*dx^3,3)./dx^2/100; % in number of molecules per voxel of 1um^2 and 100um height
%     Ctime(i,:) = interp1(dst,Cdisc(ix)',dgrid); 
% end
% Cmean = nanmean(Ctime); 
% Cmx = cummax(Ctime); 

f_Ccylinder = @(r,t) N/h/4/pi/D./t.*exp(-r.^2/4/D./t); 
[Tmesh,Dmesh] = meshgrid(Tgrid,dgrid); 
Ccylinder = f_Ccylinder(Dmesh,Tmesh); 
Ccylinder=Ccylinder'; 
Ctime = Ccylinder;
Cmx=cummax(Ctime); 

%%
figure(4)
clf

C50 = 2.9/5; % in NPE-ATP units, need to convert later... 
p=0.1:0.1:0.9; 
for i=1:9
    subplot(3,3,i)
    cm = C50*p(i)./(1-p(i)); 
    cm=cm*600; 
    contour(Tmesh,Dmesh,Ctime',[cm cm]);
    hold on
    %%
    g=0:150:1500;
    Tfirst = nan(size(g));
    for j=1:size(P,2),
        tmp = find(P(:,j)>p(i),1);
        if ~isempty(tmp),
            Tfirst(j)=Tgrid(tmp);
        end,
    end
    plot(Tfirst,g,'r')
    title(sprintf('Response fraction: %g',p(i)))
end


%%
figure(5)
clf
C50=0.29; 
A=Cmx./(Cmx+C50*600);
rmx = @(p) 2*((N*exp(-1))./(4*(C50*p./(1-p)*600)*pi*h)).^(1/2);
tmx = @(p) (N*exp(-1))./(4*(C50*p./(1-p)*600)*D*pi*h); 


% subplot(1,2,1)
plot(rmx(0:0.1:1),0:0.1:1)
ylabel('% cell responding')
xlabel('distance')

Pcmmx = cummax(P); 
Tgrid2 = linspace(0,300,251); 

for Tcomp=1000
    % Tcomp=50;
    Pexp = Pcmmx(find(Tgrid2>=Tcomp,1),:);
    Pexp(1)=NaN; % remove wound point
    clf
    hold on
    plot(0:150:1500,Pexp,'r')
    xlim([0 1500])
    plot(dgrid,A(find(Tgrid>Tcomp,1),:),'m')
    pause(0.5)
end
%
% Tmxexp=nan(11,1); 
% Tmxtheory=nan(numel(dgrid)); 
% for i=1:11,
%     [~,mi]=max(P(:,i)); 
%     Tmxexp(i)=Tgrid2(mi); 
% end
% for i=1:numel(dgrid)
%     [~,mi]=max(A(:,i));
%     Tmxtheory(i)=Tgrid(mi);
% end
% Tmxexp(end)=NaN;
% subplot(1,2,2)
% hold on
% plot(Tmxexp,0:150:1500,'r')
% plot(Tmxtheory,dgrid,'b')
% ropt = @(C50,p) 2*((N*exp(-1))./(4*(C50*p./(1-p)*600)*pi*h)).^(1/2); 
% x = lsqcurvefit(ropt,0.5,Pexp(2:end),(150:150:1500)'); 


%%
figure(1)
clf
subplot(2,2,1)
imagesc(dgrid,Tgrid,log10(Ctime/600),[-2 2])
xlabel('Distance [um]')
ylabel('Time [sec]')
ylim([0 300])
colorbar

subplot(2,2,2)
imagesc(dgrid,Tgrid,log10(Cmx/600),[-2 2])
xlabel('Distance [um]')
ylabel('Time [sec]')
ylim([0 300])
colorbar

subplot(2,2,3)
imagesc(dgrid,Tgrid,A,[0 1])
xlabel('Distance [um]')
ylabel('Time [sec]')
ylim([0 300])
colorbar

subplot(2,2,4)
imagesc(0:150:1500,Tgrid2,Pcmmx,[0 1])
colorbar



% %% 
% Cmin=10; 
% r = @(t) 4*D*t .* log(N./((4*pi*t).^(3/2))/Cmin); 
% r2 = @(t) 4*D*t; 
% Tgrid = (0:900)'; 
% plot(Tgrid,r(Tgrid),'b',Tgrid,r2(Tgrid),'r')
% figure, plot(Tgrid,r(Tgrid)./r2(Tgrid),'g')

% %% for 3D case
% syms Cmin N t r D
% r = sqrt(4*D*t*log(N/(4*pi*D*t)^(3/2)/Cmin)); 
% drdt = diff(r,t); 
% tmx = solve(drdt,t); 
% rmx = subs(r,tmx); 
% 
%% for 2D (thin cylinder) case
syms Cmin N t r D h
r=sqrt(4*D*t*log(N/h/4/D/t/pi/Cmin));
drdt = diff(r,t);
tmx = solve(drdt,t);
pretty(simplify(tmx));
rmx = subs(r,t,tmx);
pretty(simplify(rmx))

%%
