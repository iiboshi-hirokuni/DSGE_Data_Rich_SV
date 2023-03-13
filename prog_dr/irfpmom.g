/* filename:    irfpmom.g
** description: Compute posterior mean and std and CI based on
**              output of posterior simulator
*/

new;
closeall;
library user, pgraph;
dlibrary -a PszTgSen,PsdTgSen;
cls;
outwidth 128;

/******************************************************************
**         Load Parameter Draws from MH Output
*/
#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g 

mhrun  = "3";    @ 3: Case A, 4: Case B, 5: Case C   @

lpath  = resupath $+ "\\mhrun" $+ mhrun;
opath  = resupath $+ "\\mhrun" $+ mhrun;

ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "ir";
odraws = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "irs";

open fhdraw = ^ldraw for read;

nburn       = 100;        /* Number of initial draws to be discarded */
hpdprob     = 0.90 ;

drawrow = readr( fhdraw, nburn); 
drawdim = cols(drawrow);
drawrow = seekr( fhdraw, nburn);

/* Part 1: Compute the mean of x(i) and x^2(i)
*/
drawmean   = zeros(1,drawdim);
drawsqmean = zeros(1,drawdim);

eofloop = 0;
nblock  = 100;   @100@
ncols   = 20;  
ndraws  = 0;

do until eofloop; 

   drawblock = readr( fhdraw, nblock );

   drawmean  = drawmean    + sumc(drawblock)';
   drawsqmean = drawsqmean + sumc(drawblock^2)';

   ndraws = ndraws + rows(drawblock);

   locate 1,1;
   "Part 1";
   "Draws " ndraws;

   eofloop = eof(fhdraw);

endo;

drawmean   = drawmean/ndraws;
drawsqmean = drawsqmean/ndraws;
drawstdd   = sqrt(drawsqmean - (drawmean)^2);


/* Part 2: HPD Interval
*/
cls;
drawrow     = seekr( fhdraw, nburn); 
drawci      = zeros(2,drawdim);

j = 1;
do until j > drawdim/ncols;

   /* Read only the j'th set of columns
   */
   drawrow   = seekr( fhdraw, nburn); 
   drawblock = readr( fhdraw, 1);

   drawcol = drawblock[1,(j-1)*ncols+1:j*ncols];
   eofloop = 0;
   ndraws  = 1;
   
   do until eofloop; 

      drawblock = readr( fhdraw, nblock );
      drawcol   = drawcol | drawblock[.,(j-1)*ncols+1:j*ncols];
      ndraws    = ndraws + rows(drawblock);
      eofloop   = eof(fhdraw);

      locate 1,1;
      "Part 2";
      "ndraws" ndraws;
      "Columns" (j-1)*ncols+1 "to" j*ncols;
         
   endo;

   l = 1;
   do until l > ncols;  
 
      locate 4,1;
      "HPD Intv for Column" (j-1)*ncols+l;
      drawci[.,(j-1)*ncols+l] = hpdint(drawcol[.,l],hpdprob);

      l = l+1;
   endo;

   j = j+1;

endo;
clear drawcol;

create fhdraws = ^odraws with DRAWSUM, drawdim, 8;
writer(fhdraws, real(drawmean | drawstdd | drawci) );
closeall fhdraws; 

end;

/**************************************************************/
/*                         Procedures                         */
/**************************************************************/

proc(1) = hpdint(draws,percent);
local     ndraws, drawdim, hpdband, nwidth, i, drawcoli,
          bup, minwidth, newwidth, j;

drawdim   = cols(draws);
ndraws    = rows(draws);
hpdband   = zeros(2,drawdim);
nwidth    = int(percent*ndraws);

i = 1;
do until i > drawdim;
   drawcoli = draws[.,i];
   /* sort response for period i, element 1 is max
   */
   drawcoli = rev(sortc(drawcoli,1));
   bup   = 1;
   minwidth  = drawcoli[1] - drawcoli[nwidth];
   j = 2;
   do until j > (ndraws-nwidth+1);
      newwidth = drawcoli[j] - drawcoli[j+nwidth-1];
      if newwidth < minwidth;
         bup = j;
         minwidth = newwidth;
      endif;
      j = j+1;
   endo;
   hpdband[2,i] = drawcoli[bup];
   hpdband[1,i] = drawcoli[bup+nwidth-1];
   i = i+1;
endo;
retp(hpdband);
endp;
