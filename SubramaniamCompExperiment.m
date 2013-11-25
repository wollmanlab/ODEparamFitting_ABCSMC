%% Computational experiments on the Subramaniam Model  
%% Run original model
ver = 0;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% Run version 1
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% kf3 at 10%
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
K(38) = K(38)/10;
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% kf3 at 1%
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
K(38) = (K(38)/100);
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% Km_cai_3 x10
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
K(39) = K(39)*10;

SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% kf4 at 10%
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
K(40) = K(40)*10;
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% Receptor at 10 time original
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
x0(2) = x0(2)*10;
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% kf3 tune down to 1%, Kmcai3 x 100, kf4 down to 1%
% kf3 manipulation
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
K(38) = K(38)/10;
K(39) = K(39)*10;
K(40) = K(40)/10;
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
%% Turn up receptor  to x 100
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
x0(2) = x0(2)*100;
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);
SubramaniamConservation(K,x0,L,t0,t,dt);
%% Turn up the ip3 generation reaction in flux 18
ver = 1;
[K,x0,L,t0,t,dt] = SubramaniamInput(ver);
K(57) = K(57)*10;
SubramaniamModelTimeCourse( K,x0,L,t0,t,dt);