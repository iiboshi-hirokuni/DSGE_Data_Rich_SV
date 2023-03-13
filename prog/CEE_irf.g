/*******************************************************/
/*                                                     */
/* Computing Impulse responses based on                */
/* parameters of DSGE models                           */
/*                                                     */
/*******************************************************/

/* filename:    nkbc_irf.g
** description: The program converts the DSGE parameters
**              into impulse response functions
** created:     21/10/99
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

#include c:\dge_sw\pathspec.g;    @@
#include c:\dge_sw\ceespec.g;    @@


nirf = 20;    @number of Horizon@

/* Open files with Parameter Draws
*/
mhrun  = "3";       @  2: SW model, 3: Case A, 4: Case B, 5: Case C   @

lpath  = parapath $+ "\\mhrun" $+ mhrun;
opath  = resupath $+ "\\mhrun" $+ mhrun;

lpara  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa";
open fhpara = ^lpara for read;

/* Define some Constants
*/
nblock  = 10000;    /* block size */
nuse    = 5;    /* use every nuse observation */

/* Initialize Output files
*/
oirf = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "ir";
create fhirf = ^oirf with IRF, nvar*nshock*nirf, 8;

   
/* Initialization
*/ 
eofloop = 0;
loopct  = 1;

do until eofloop; 

   /* Read blocks of size nblock
   */
   parasim = readr( fhpara, nblock );
   nread = rows(parasim);

   indseq = seqa(1,1,nread);
   indseq = indseq % nuse;
   indseq = indseq .NE 0;
 
   parasim = delif(parasim,indseq);
   nread = rows(parasim);

   irfsim = zeros(nread,nvar*nshock*nirf);
   
   j = 1;   
   do until j > nread;
   
      /* For each element of the block, compute the desired statistic
      */
      para = parasim[j,1:27]';

      ssirf = dsgeirf(para,nirf);
     
      irfsim[j,.] = vec(ssirf)';
      
      j = j+1;
      
   endo;

   /* Write Output
   */


 locate 1,1;
   "Block" loopct;
   "Draws Processed:" nread;

     wr = writer(fhirf,irfsim);

   loopct = loopct + 1;
   eofloop = eof(fhpara);


endo;

closeall fhpara, fhirf;

/**********************************************/

proc (1) = dsgeirf(para,nirf);
/* INPUT:
**   para   : parameter vector 
**   nirf   : length of irf
** OUTPUT:
**   ssirf  : impulse response 
**            the structural shocks.
**            0 if retcode=0
**            (nirf by 25)
*/
local sig_cov, sig_chol, impact,
      T1,TC,T0,TETA,RC,RR, ZZ,
      yyirfall, yyirf, ss, sh_ind, t;

/* Compute Covariance Matrix 
** and its Chol decomposition
*/ 
sig_chol = eye(9);        @@
sig_chol[9,9] = 1/4;             @ monetary    @

/*
sig_chol[1,1] = para[19];             @ preference @
sig_chol[2,2] = para[20];             @ investment @
sig_chol[3,3] = para[21];             @  Tobin's q @
sig_chol[4,4] = para[22];             @ labor      @
sig_chol[5,5] = para[23];             @ wage markup @
sig_chol[6,6] = para[24];             @ productivity      @
sig_chol[7,7] = para[25];             @ price markup @
sig_chol[8,8] = para[26];             @ goverment   @
sig_chol[9,9] = para[27];             @ monetary    @
*/

/* solve DSGE
*/
{T1,TC,T0,RC} = dsgesolv(para);

yyirfall = zeros(nirf,nvar*nshock); 

sh_ind = 1;
do until sh_ind > nshock;
   
   impact = sig_chol[.,sh_ind];
   yyirf  = zeros(nirf,nvar);

  ZZ      = zeros(9, rows(T1));                          @@
ZZ[1,1] = 1;               /* output gap */
ZZ[2,7] = 1;               /*  consumption */
ZZ[3,2] = 4;               /* inflation */
ZZ[4,8] = 4;               /*  nominal rate  */
ZZ[5,3] = 1;               /* wage */
ZZ[6,9] = 1;              /* cost of capital */
ZZ[7,4] = 1;               /* capital stock */ 
ZZ[8,6] = 1;               /* investment */            
ZZ[9,10] = 1;              /*  labor */   

   ss     = (T0*impact);
   yyirf[1,.] = (ZZ*ss)';
      
   t = 2;
   do until t > nirf;
      ss = T1*ss;     
      yyirf[t,.] = (ZZ*ss)';
      t = t+1;
   endo;

  @ yyirf[.,1] = cumsumc(yyirf[.,1]); @  /* output response is cumulative */

     yyirfall[.,1+nvar*(sh_ind-1):nvar*sh_ind] = yyirf;

   sh_ind = sh_ind+1;
   
endo;

retp(yyirfall);
endp;
      
#include c:\dge_sw\ceemod.src;    @@

