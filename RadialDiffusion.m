%{ 
this function implements the finite difference methods to simulate  radially symmetric diffusion
inputs:
       D = diffusion coefficient in micrometer^2/sec
       IC = the initial condition of the spatial vector, which is a n x1
       vector
       dr = the mesh spacing in micrometer
       dt = the time step, in seconds
       tNum = the duration of simulation, in number of time steps (not including the initial condition)
       tridiagFlag = a flag to indicate whether to use the tridiag function
output:
       m = a matrix whose rows represent the time and the columns represent
       the spatial grid points. The entries represent the concentration at
       a specific grid point at a speciic time


%}
function m = RadialDiffusion(D,IC,dr,dt,tNum,tridiagFlag)
   % calculate the B constant to be used in the  
   B = D*dt/(2*dr^2);
   m = zeros(length(IC),tNum+1);
   m(:,1) = IC;
   %% implicit
   if tridiagFlag == 1  
       % prepare the inputs for the tridiag function
       a = zeros(length(IC),1);
       a(1) = -1/dr;
       a(end) = 1;
       a(2:end-1) = 1+2*B;
       b = zeros(length(IC),1);
       b(1:end-1) = -B;
       b(end) = 0;
       c = zeros(length(IC)-1,1);
       c(1) = 1/dr;
       c(2:end) = -B;
       for i = 1:tNum
           % This is a tweak to set the two concentration equal so that no
           % molecules flow to the left of the second grid point
           %if i > 2
            %   m(1,i) = m(2,i);
           %end
           % compute the f vector for the tridiag function
           disp(strcat('tNum = ',num2str(i)));
           f = zeros(length(IC),1);
           f(1) = 0; % the left boundary has zero flux
           f(end) = 0; % the right boundary has zero concentration
           parfor ind = 2:length(f)-1
               f(ind) = (B-D/(2*ind*dr))*m(ind-1,i) + (1-2*B)*m(ind,i) + (B + D/(2*ind*dr))*m(ind+1,i); 
           end
           m(:,i+1) = tridiag(a,b,c,f);
       end
   % explicit method implementation
   else
      A = sparse(length(IC),length(IC));
      % code the elements of the sparse matrix 
      for i = 2:length(IC)-1
          B = D*dt/(dr^2);
          C = 1/(2*i*dr);
          A(i,i-1) = (B -C);
          A(i,i) = (1 - 2*B);
          A(i,i+1) = (B + C);
      end
      % iterate through to calculate the time evolution 
      for i = 1:tNum
          disp(strcat('tNum = ',num2str(i)));
          m(:,i+1) = A*m(:,i);
      end
       
   end
end