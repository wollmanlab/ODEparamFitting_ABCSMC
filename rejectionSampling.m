function [FinalPop,FinalScores,FinalSim,TotalAttempts] = rejectionSampling(fScore,fRand,fEval,TargetSize,varargin)

arg.verbose = false; 
arg = parseVarargin(varargin,arg); 
t0=now; 

FinalScores = [];
FinalPop=[]; 
FinalSim={}; 
counter=0; 
while numel(FinalScores)<TargetSize 
    counter=counter+1; 
    if counter==1 || isempty(FinalScores)
        TotalAttempts(counter)=TargetSize*4; % assume 25% success rate is as good as it gets; 
    else
        TotalAttempts(counter)=min(960,ceil(sum(TotalAttempts)*TargetSize/numel(FinalScores)/4)); % quarter of what we need to accomplish
    end
    Sim=cell(TotalAttempts(counter),1);    
    Scores=inf(TotalAttempts(counter),1);
    Pop=fRand(TotalAttempts(counter));
%     fprintf('inner iter %g function calls: %g',counter,TotalAttempts(counter))
    parfor i=1:TotalAttempts(counter)
        warning off
        [Scores(i),Sim{i}] = fScore(Pop(i,:));
    end
    ix = fEval(Scores,Pop); 
    FinalPop=[FinalPop; Pop(ix,:)];  %#ok<*AGROW>
    FinalScores=[FinalScores; Scores(ix)]; 
    FinalSim=[FinalSim; Sim(ix)]; 
    arg.verbose && fprintf(' total function calls: %g. overall time %s, pop size is %g pop mean is: %g\n',sum(TotalAttempts),datestr(now-t0,13),numel(FinalScores),nanmean(FinalScores));
end
   
FinalPop=FinalPop(1:TargetSize,:); 
FinalScores=FinalScores(1:TargetSize); 