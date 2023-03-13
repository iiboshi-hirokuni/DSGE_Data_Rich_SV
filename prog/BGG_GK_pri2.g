/* filename:    BGG_GK_prior2.g
** description: The program deletes the draws that lie outside of the parameter space
** created:     12/07/2010
*/

new;
closeall;
library pgraph, user;
dlibrary -a PszTgSen,PsdTgSen;

cls;

/* Global defaults 
*/
  __output = 1;
  _fcmptol = 1E-10;

#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g

/* Define some Constants
*/
nsim   = 100;   /* block size */
mhrun  = "0";

/* Open files with Parameter Draws
*/
lpath  = parapath $+ "\\mhrun" $+ mhrun;  

ltheta = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa0";
open fhltheta = ^ltheta for read;

/* Initialize Output files
*/
opath = lpath;

otheta = opath  $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa";
create fhotheta = ^otheta with THETA, npara, 8;

odet  = opath  $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "dt";
create fhdet = ^odet with DET, 2, 2;

/* Initialization
*/ 
eofloop = 0;
loopct  = 0;
nread   = 0;
nvalid  = 0;

do until eofloop; 

   /* Read blocks of size nsim
   */
   thetasim = readr( fhltheta, nsim );  

   nread   = nread + rows(thetasim);

   /* Eliminate draws that are in the non-existence and indeterminacy region
   */
   invaliddraw = zeros(rows(thetasim),1);
   detsim      = zeros(rows(thetasim),2);

   j = 1;
   do until j >  rows(thetasim);

      /* For each element of the block, compute the desired statistic
      */

      {T1,TC,T0,RC} = dsgesolv(thetasim[j,.]');
      
      if (RC[1] == 1) AND (RC[2]==0);
         /* indeterminacy */
         detsim[j,2] = 1;
         
      elseif (RC[1] == 1) AND (RC[2]==1); 
         /* determinacy */
         detsim[j,1] = 1;   
       
      else;
         invaliddraw[j,1] = 1;  

      endif; 

      j = j+1;
  endo;

  RC; 

  thetasim = delif(thetasim, invaliddraw);  
  detsim   = delif(detsim, invaliddraw);     //cols(thetasim); invaliddraw;
  thetasim = delif(thetasim, detsim[.,2]);

  
   nvalid  = nvalid + rows(thetasim);

  /* Write Output
  */
  wr = writer(fhdet,detsim); 
  wr = writer(fhotheta,thetasim);
  
  loopct = loopct + 1;

  locate 1,1;
  "Loop            " loopct;
  "Relative Size of Regions [det, indet]"; meanc(detsim)';
  eofloop = eof(fhltheta);
  
endo;

"Fraction of valid draws" nvalid/nread;

closeall fhdet, fhltheta, fhotheta;  


#include c:\dge_BGG_GK\BGG_GK_mod.src

