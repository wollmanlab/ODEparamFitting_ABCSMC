%% This function genertes the input parameters for the spatiotemporal simulation

function  [dr,xDim,yDim,dt,tNum,totalLig,deg,reaction,D] = SpatioTemporalInput(woundSize,cellDens)
   %{
       The default units for all variables are in micrometer and seconds (the case when spatial conversion is 1 pixel = 1 micrometer and time conversion is 1 step = 1 sec )
       input:
             woundSize = the diameter of the wound in micrometer
             cellDens = cell density. The # of cells in a 2000*2000 micron
             squared (which is the usual size of the image)

       output:
         woundSize = diameter of the wound
         height = height of the chamber
         D = diffusion coefficient of ATP
         numATP = total number of ATP storage per cell
         h = height of cylinder 
         numCell = number of cells in the simulation
         xDim = x dimension of the simulation 
         yDim = y dimension of the simulation 
         dt = size of time step
         dDim = unit of spatial increment
         IC = initial condition vector  
         mode = mode of ATP release. 1 = dilution, 2 = diffusion with
         degradation, 3 = wave propagation with positive feedback
         r is the radius of the simulation
         xDim = x dimension in units of spatial steps
         yDim = y dimension in units of spatial steps
         deg = degradation rate
         reaction = a function that describes the ligand incuded ligand
         secretino by the cells. This is used in a positive feedback loop. It is modeled by a form of Hill equation 
     %}
   % the spatial conversion variable is in unit of pixel/micrometer
   spatialConv = 0.5;
   % the time conversion variable is in unit of time step/sec
   timeConv = 10;
   % cellDens is the spatial density of the cell in units of cells/ micrometer squared 
   cellDens = cellDens/(2000^2)/(spatialConv^2);
   woundSize = woundSize*spatialConv;
   D = 200*(spatialConv)^2/timeConv;
   ligPerCell = 1e14;
   totalLig = (pi*(woundSize/2)^2*cellDens*ligPerCell); 
   dt = 1/timeConv;
   tNum = 600/dt;
   dr = 1/spatialConv;
   % ligand hydrolysis constant. The units are in # ligand*micronsquared/(cell*sec). It is directly proportional to cell density. 
   deg = 1*spatialConv^2/timeConv;
   % maximum rate of ligand production. the units are in # ligands*(micron
   % squared)/(cell*sec)
   Vmax = 1e9*(spatialConv)^2/timeConv;
   %dimensions of the simulation 
   xDim = 2000*spatialConv;
   yDim = 2000*spatialConv;
   % n = Hill coefficient
   n = 2;
   % Km = equilibrium constant in units of ligand/(micrometer squared)
   Km = 10e2/(spatialConv^2);
   reaction = @(ligand) Vmax*cellDens*(ligand.^n/(Km.^n + ligand.^n));
   %reaction = @(ligand) Vmax*cellDens*ligand;
end