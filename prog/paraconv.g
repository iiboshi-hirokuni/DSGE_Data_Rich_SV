/* filename:    paraconv.g
** description: Informal assessment of convergence:
**              plot recursive and moving averages
*/

new;
closeall;
library user, pgraph;
cls;
outwidth 128;

mhrun  = "0";       @ 0: prior   2: posterior @

nburn       = 1;        /* Number of initial draws to be discarded */
nblock      = 100;
nma         = 10;
nstart      = 10;

/******************************************************************
**         Load Parameter Draws from MH Output
*/
#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g

lpath  = parapath $+ "\\mhrun" $+ mhrun;
ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa";

//load path = ^priopath para_names;

para_names={kappa h sigma_c sigma_L phi iota_p iota_w  
             theta_p theta_w rho_r  mu_pi  mu_y  
             rho_a   rho_c   rho_k  rho_e  rho_f  rho_g  rho_l
             e_a   e_c  e_e  e_f  e_g e_k e_l e_r };

open fhdraw = ^ldraw for read;

drawrow = readr( fhdraw, nburn); 
drawdim = 27;   //cols(drawrow);
drawrow = seekr( fhdraw, nburn);

fhdraw;


graphset;
fonts("microb");
_ptitlht = 0.4;
_paxht = 0.4;
_pnumht = 0.4;
_pltype = {6 1 2}; 
_plegstr = "Cumulative Average\000Moving Average"; 
_plegctl = {2 4 4 4};
_plwidth = 3.5;

begwind;
window(5,6,0);    @@

j = 1;
do until j > drawdim;

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
      eoflopp = 1;
               
   endo;
   
   draw_cummean = cumsumc(drawcol)./seqa(1,1,ndraws);
   draw_ma      = movave(drawcol,nma);
   draw_ma      = draw_cummean[1:nma-1,. ] | draw_ma;
   draw_x       = seqa(1,1,ndraws);

  title(para_names[j]);
   xy(draw_x[nstart:ndraws,.],draw_cummean[nstart:ndraws,.]~draw_ma[nstart:ndraws,.]);
   nextwind;

   "Generated Graph" j "out of" drawdim;

   j = j+1;
   clear drawcol;
   
endo;

endwind;

closeall fhdraw;

end;

/**************************************************************/
/*                         Procedures                         */
/**************************************************************/

proc(1) = movave(draws,nma);
local drawdim, ndraws, mavec, i, j;    

drawdim   = cols(draws);
ndraws    = rows(draws);

mavec = zeros(1+ndraws-nma,drawdim);

i = 1;
do until i > drawdim;

   j = 1;
   do until j > (1+ndraws-nma);
      mavec[j,i] = meanc(draws[j:j+nma-1,i]);
      j = j+1;
   endo;

   i = i+1;
endo;

retp(mavec);
endp;
