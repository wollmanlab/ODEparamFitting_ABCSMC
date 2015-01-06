function status = timelimit(~,~,flag)
maxRunTime=5; % seconds
persistent starttime; 
if strcmp(flag,'init')
    starttime=cputime; 
end
runtime=cputime-starttime; 

if runtime>maxRunTime
    status=1;
%     fprintf('\nhalting execution of ODE, runtime: %g\n',runtime)
else
    status=0; 
end