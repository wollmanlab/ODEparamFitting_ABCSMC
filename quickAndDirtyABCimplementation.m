%% Init
clear 
close all
clc

%% include
addpath /home/jason/calciummodeling/
addpath /home/jason/calciummodeling/Traj_Compare_Conditions/Traj_Orig_Ca/

%% load experimental data
load /home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/traj_147.mat
load /home/jason/calciummodeling/Traj_Folder_L_10_Run1_Feb232014/tInterp.mat
L=10; 
%% Get the initial population
input=BennettModelInput(1); 
input(17:19)=[]; 
input=input(:)';

PopSize = 25;
Eps = 150;
InitialPop = [];
InitialLogPop = [];
tNow = now;
PopScore = [];
Score = inf;
% define the new function to calculate score
fOpt = @(norminput) BennettABCOptimize(10.^(norminput).*input,L,tInterp,traj_interp);
i = 0;
while size(InitialPop,1)<PopSize
    % the flag is set for valid simulation. The flag is 1 if valid
    flag =0;
    logInput = zeros(1,length(input));
    while flag==0
        %disp('new input');
        logInput = rand(1,32) -0.5*4;
        logInput(7) = rand(1,1);
        [Score,caSim,g] = fOpt(logInput);
        % check equality constaints
        if Score == inf
            flag =0;
        else
            flag =1;
        end
        %if any(g(1:14)) || g(15) > 100 || g(16) < 0.03 || g(16) > 0.1 || g(17)< 0.5 || g(17)>0.8
        %    flag =0;
        %else 
        %    flag =1;
        %end
    end
    if Score < Eps
        InitialLogPop = [InitialLogPop;logInput];
        InitialPop = [InitialPop; 10.^(logInput).*input];
        PopScore = [PopScore; Score];
        fprintf('iter: %g, time: %s, popsize: %g\n',i,datestr(now-tNow,13),size(InitialPop,1));
    end
end
CurrentGoodPop = InitialLogPop;
% save the initial population and the population score
save('InitialPop.mat','InitialPop');
save('InitialLogPop.mat','InitialLogPop');
save('InitialScore.mat','PopScore');
%% reapeat for a series of eps schedule
EpsVec=[100 33.3 10 3.33 1 0.33 0.1];
Score = inf;
for i=1:numel(EpsVec)
    Eps=EpsVec(i);
    NewLogPop = zeros(size(InitialLogPop));
    NewPop=zeros(size(InitialLogPop));
    PopScore=zeros(size(InitialLogPop,1));
    
    % Perturb each one of the initial population
    parfor j =1:PopSize
        flag = 0;
        Score = inf;
        newMember = zeros(1,length(input));
        
        while flag == 0 || Score >= Eps;
            % perturb the current member of the population
            newMember = CurrentGoodPop(j,:) + randn(1,32)/2;
            newMember(7) = min(newMember(7),1);
            newMember(7) = max(newMember(7),0);
            [Score,caSim,g] = fOpt(newMember);
            if Score == inf
                flag =0;
            else
                flag =1;
            end
            % check equality constaints
            %if any(g(1:14)) || g(15) > 100 || g(16) < 0.03 || g(16) > 0.1 || g(17)< 0.5 || g(17)>0.8
            %    flag =0;
            %else 
            %    flag =1;
            %end
            
        end
        NewLogPop(j,:) = newMember;
        NewPop(j,:) = 10.^(newMember).*input;
        PopScore(j) = Score;
        fprintf('epsIter: %g, time: %s, popsize: %g\n',i,datestr(now-tNow,13),size(NewPop,1));
        disp(strcat('j = ',num2str(j)));
            
    end
    % save current population
    save(strcat('Pop',num2str(i)),'NewPop');
    save(strcat('LogPop',num2str(i)),'NewLogPop');
    save(strcat('PopScore',num2str(i)),'PopScore');
    
    CurrentGoodPop = NewLogPop;
    fprintf('Finished an outter iteration: %g time: %s - yay!\n',i,datestr(now-tNow,13));  
end


