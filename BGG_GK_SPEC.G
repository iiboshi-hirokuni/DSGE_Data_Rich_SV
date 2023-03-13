/* Model Specification
**
** 1 = central bank reacts to inflation, with interest rate smoothing
** 2 = central bank reacts to inflation, no interest rate smoothing
*/
mspec  = 1;
lmodel = "m" $+ ftos(mspec,"%*.*lf",1,0) ;

mprior = 0;
lprior = ftos(mprior,"%*.*lf",1,0);

/*
** para_mask is npara * (number of model specifications)
**   to fix a parameter under specification mspec, set the corresponding
**   element of para_mask to 1
** para_fix is npara * (number of model specifications)
**   and contains the values for the fixed parameters (all other entries
**   are irrelevant
*/
nvar   = 11;  nshock = 8;
npara  = 38;
ngrid  = 100;  /* number of grid points for density plots */


             
/* Selected data set
** 0 = actual, full sample
** i = simulated data set i
*/
datasel    = 2;    @2:USA, 1:Simulation @
 dataselstr = ftos(datasel,"%*.*lf",1,0);


if datasel ne 1; 

   subT = 1;   
   
   subTstr = ftos(subT,"%*.*lf",1,0);
   
elseif datasel == 1; /* simulated data */

   Tsimtotal = 200;
   Tsimburn  = 200;
   Tsim = 200 | 500;
   subT = 1;   
   subTstr = ftos(subT,"%*.*lf",1,0);
     
endif;
