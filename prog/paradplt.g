/* filename:    paradplt.g
** description: Plot prior and posterior densities:
**              
*/

new;
closeall;
library user, pgraph;
cls;
outwidth 128;


/******************************************************************
**         Load Density Estimates
*/
#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g

npara = 27;

mhrun   = "0";
lpath   = parapath $+ "\\mhrun" $+ mhrun;
ldens0  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pad";


mhrun   = "0";
lpath   = parapath $+ "\\mhrun" $+ mhrun;
ldens1  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pad";


open fhdens0 = ^ldens0 for read;
open fhdens1 = ^ldens1 for read;

densmat0 = readr( fhdens0, ngrid); 
densmat1 = readr( fhdens1, ngrid); 

/* Load parameter names
*/
 @load path = ^priopath para_names;@

 para_names = "kappa" $| "h" $| "sigma_c" $| "sigma_L" $|  "Phi "$| "Iota_P"$| "Iota_W"$|
              "Theta_P" $|  "Theta_W" $| "Rho_R" $| "Mu_Pi" $| 
              "Mu_Y" $| "Rho_A" $|  "Rho_C" $|  "Rho_K" $| "Rho_E" $|  "Rho_F" $| 
              "Rho_G" $| "Rho_L" $|
              "e_A" $|  "e_C" $|  "e_E" $|  "e_F" $|  "e_G" $|  "e_K" $|  "e_L" $|  "e_R"   ;  

 

/* Load gridrange for density evaluation
*/

gridfile = "paragrid.out";
load path=^priopath gridrange[npara,2] = ^gridfile;



l =1;

do until l >27;
j = 1;
graphset;
fonts("microb");
_ptitlht = 0.4;
_paxht = 0.4;
_pnumht = 0.4;
_pltype = {6 1 2}; 
_plegstr = "Prior\000Posterior"; 
_plegctl = {0 4 4 4};
_plwidth = 3.5;

begwind;
window(4,3,0);      @@

do until (j > 12)or(l>27);

   /* Compute Grid for density evaluation */  rows(gridrange);
   paragrid = seqa(gridrange[l,1],(gridrange[l,2]-gridrange[l,1])/(ngrid-1),ngrid);

  
    title(para_names[l]);
     
   xy(paragrid,densmat0[.,l]~densmat1[.,l]);
   @  xy(paragrid,densmat0[.,l]); @

   nextwind;

   "Generated Graph" l "out of" npara;

   j = j+1;
   l = l + 1;
   
   clear drawcol;
   
endo;

endwind;  
  

"Press Key to Continue";
k = keyw;

  //k = K +1;
endo;




/*

 dens= zeros(ngrid*3,2);

j = 1;
do until j > 3;

   dens[1+ngrid*(j-1):ngrid*j,1] = densmat0[.,j];

j= j+1;
endo

OUTPUT FILE=c:\cee-dge\dens.OUT RESET;

 densmat0;

 OUTPUT OFF;

*/

closeall fhdens0, fhdens1;

end;

