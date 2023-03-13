/* filename:    paradens.g
** description: Computes posterior density estimates 
**              based on the output of the posterior simulator
*/

new;
closeall;
library user, pgraph;
cls;
outwidth 128;

mhrun  = "0";       @ 0: prior,    2: posterior   @

nburn   = 2;        /* Number of initial draws to be discarded */
nblock  = 10;

/******************************************************************
**         Load Parameter Draws from MH Output
*/
#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g

npara = 27;

lpath  = parapath $+ "\\mhrun" $+ mhrun;
opath  = parapath $+ "\\mhrun" $+ mhrun;

ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa";
odens  = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pad";

open fhdraw = ^ldraw for read;

drawrow = readr( fhdraw, nburn); 
drawdim = cols(drawrow);
drawrow = seekr( fhdraw, nburn);

/* Load gridrange for density evaluation
*/
gridfile = "paragrid.out";
load path=^priopath gridrange[npara,2] = ^gridfile;

/* initialize matrix with density
** estimates
*/
densmat = zeros(ngrid,npara);

j = 1;
do until j > drawdim;

   locate 1,1;
   "Processing column" j;

  @   if para_maskinv[j,mspec] == 1; @

      /* Parameter has not been fixed
      */

      /* Read only the j'th column
      */
      drawrow   = seekr( fhdraw, nburn); 
      drawblock = readr( fhdraw, 1);

      drawcol = drawblock[1,j];
      eofloop = 0;
      ndraws  = 1;
   
      do until eofloop; 

         drawblock = readr( fhdraw, nblock );
         drawcol   = drawcol | drawblock[.,j];
         ndraws    = ndraws + rows(drawblock);
         eofloop   = eof(fhdraw);
         
      endo;

      /* Compute SD for bandwidth selection */
      drawstdd = stdc(drawcol);
      
      /* Compute Grid for density evaluation */
      paragrid = seqa(gridrange[j,1],(gridrange[j,2]-gridrange[j,1])/(ngrid-1),ngrid);
      
      /* Construct density estimate */
      h    = 1.06*ndraws^(-1/5)*drawstdd;
      densmat[.,j] = densest(paragrid,drawcol,h);      
      
    @ endif;   @
     
   j = j+1;

endo;
clear drawcol;

/* Save the results
*/
cls;

create fhdens = ^odens with PARA, npara, 8;
wr = writer(fhdens, densmat );
closeall fhdraw, fhdens; 

"Saved Density Estimates";

end;

/**************************************************************/
/*                         Procedures                         */
/**************************************************************/


proc (1) = densest(paragrid,draws,h);
/* This procedure computes a simple kernel density estimate
*/
local   ndraws, npt, pdfhat, i;

ndraws = rows(draws);
npt    = rows(paragrid);
pdfhat = zeros(npt,1); 

i = 1;
do until i > npt;
   pdfhat[i] = 1/(ndraws*h)*sumc( pdfn( (ones(ndraws,1)*paragrid[i]-draws)/h ) );
   i = i+1;
endo;

retp(pdfhat);
endp;