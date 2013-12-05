%% This function genertes the input parameters for the spatiotemporal simulation
function  [s,height,D,numATP,numCell,r,dr,tNum,dt,mode] = SpatioTemporalInput(ver)
   %{
       
       output:
         s = diameter of the wound
         height = height of the chamber
         D = diffusion coefficient of ATP
         numATP = total number of ATP storage per cell
         h = height of cylinder 
         numCell = number of cells in the simulation
         xDim = x dimension of the simulation 
         yDim = y dimension of the simulation 
         dt = size of time step
         dDim = unit of spatial increment
         mode = mode of ATP release. 1 = dilution, 2 = diffusion with
         degradation, 3 = wave propagation with positive feedback
         r is the radius of the simulation
     %}
   s = 150;
   height = 100;
   D = 20;
   numATP = 1e11; 
   numCell = 100;
   r = 20000;
   dr = 0.1;
   tNum = 6000;
   dt = 0.1;
   mode = 1;
    
end