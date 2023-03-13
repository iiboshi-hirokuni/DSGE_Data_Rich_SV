/* filename:    varipmom.g
** description: Compute posterior mean and std and CI based on
**              output of posterior simulator
*/

new;
closeall;
library user, pgraph;
cls;
outwidth 128;

#include c:\dge_sw\pathspec.g    
#include c:\dge_sw\ceespec.g   


mhrun  = "3";        @ 3:data_rich_A,  4: data_rich_B,  5: data_rich_C   @

nburn       = 10000;        /* Number of initial draws to be discarded */
hpdprob     = 0.90;       /*  Credible Band  */  
nblock      = 10000;       /*  number of draws for one block */

nobs = 60;


/******************************************************************
**         Load Parameter Draws from MH Output
*/
#include c:\dge_sw\pathspec.g;
#include c:\dge_sw\ceespec.g;

lpath  = parapath $+ "\\mhrun" $+ mhrun;
//opath  = parapath $+ "\\mhrun" $+ mhrun;

ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "sck";
//odraws = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pas";

open fhdraw = ^ldraw for read;

drawrow = readr( fhdraw, nburn); 
drawdim = cols(drawrow);
drawrow = seekr( fhdraw, nburn);

/* Part 1: Compute the mean of x(i) and x^2(i)
*/
drawmean   = zeros(1,drawdim);
drawsqmean = zeros(1,drawdim);
drawcross  = zeros(drawdim,drawdim);

eofloop = 0;
ndraws  = 0;
i = 1;

do until eofloop; 
   
   
   drawblock = readr( fhdraw, nblock );
 
   drawmean  = drawmean    + sumc(drawblock)';
   drawsqmean = drawsqmean + sumc(drawblock^2)';
   drawcross  = drawcross  + drawblock'*drawblock;  

   ndraws = ndraws + rows(drawblock);
   
   locate 1,1;
   "Part 1";
   "Draws " ndraws;

   eofloop = eof(fhdraw);

endo;

drawmean   = drawmean/ndraws;
drawsqmean = drawsqmean/ndraws;
drawstdd   = sqrt(drawsqmean - (drawmean)^2);
drawcov    = drawcross/ndraws - drawmean'*drawmean;


/* Part 2: HPD Interval
*/

@cls;@
drawrow     = seekr( fhdraw, nburn); 
drawci      = zeros(2,drawdim);
drawmean    = zeros(1,drawdim);

j = 1;
do until j > drawdim;

   // Read only the j'th column
   
   drawrow   = seekr( fhdraw, nburn); 
   drawblock = readr( fhdraw, 1);

   drawcol = drawblock[1,j];
   eofloop = 0;
   ndraws  = 1;
   i = 1;

   do until eofloop; 

      drawblock = readr( fhdraw, nblock );
  
      drawcol   = drawcol | drawblock[.,j];
      ndraws    = ndraws + rows(drawblock);
      eofloop   = eof(fhdraw);
            
         locate 1,1;
         "Part 2";
         "ndraws" ndraws;
         "Column" j; 
         
   endo;

   drawmean[.,j] = meanc(drawcol[.,1]);   
   drawci[.,j] = hpdint(drawcol,hpdprob);
   j = j+1;

endo;
clear drawcol;


/* Report Posterior Mean and stderror */
ss_mean = zeros(nobs-1,9);
ss_low  = zeros(nobs-1,9);
ss_up   = zeros(nobs-1,9);

i = 1;   
do until i > 9;
   ss_mean[.,i] = drawmean[1,1+(nobs-1)*(i-1):(nobs-1)*i]';
   ss_low[.,i] =  drawci[1,1+(nobs-1)*(i-1):(nobs-1)*i]';
   ss_up[.,i] =  drawci[2,1+(nobs-1)*(i-1):(nobs-1)*i]';

   i = i +1 ;
endo;   

cls;
"MODEL       " lmodel;
"PRIOR       " lprior;
"MH-Run      " mhrun;
"Subsample   " subT;
"Coverage    " hpdprob;  
"Dataset     " dataselstr;
"Draws " ndraws;
"=============================================";
" " ;


/****************************
 Drawing Graph of Data  
***************************/

titlestr = "Preference_Shock " $|
           "Investment_Shock " $| 
           "Equity_Premium_Shock " $|
            "Labor_Supply_Shock " $|  
           "Wage_Markup_Shock " $| 
           "Productivity_Shock " $| 
            "Price_Markup_Shock " $|                  
           "Government_Spending_Shock " $|
           "Monetary_Policy_Shock "  ; 


graphset;
    begwind;
    margin(0,0,0.2,0.2);
    window(3,3,0);           @@
    fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    _pcolor = {2 3 3 3};
    _pltype = {6 1 1 1};
    _protate = 0; 
    _plctrl = 0;
    _plwidth = 5.5;
    _psymsiz = 5;
    _plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";

       
       y_ind = 1;

       do until y_ind > 9;
       /* iterate over endogenous variables */
          
          title(titlestr[y_ind]);

          if (datasel ==2 )or (datasel == 3) ;
              xtics(1981, 1996,2,4);
              x=seqa(0/4,1/4, rows(ss_mean))+1980;
          else;
              xtics(1947, 2005,2,4);
              x=seqa(0/4,1/4, rows(ss_mean))+1947;
          endif;
          
          xy(x, ss_mean[.,y_ind]~ss_low[.,y_ind] ~ss_up[.,y_ind]);

          nextwind;  
                 
          y_ind = y_ind+1;
		  
		  if (y_ind == 5)or (y_ind == 7);
		      y_ind = y_ind+1;
		  endif;	  
		  
       endo;
       
    
    endwind;   

/* save  data */

 yy_save = {};
 y_ind = 1;
do until y_ind > 9;
   yy_save = yy_save~ss_mean[.,y_ind]~ss_low[.,y_ind] ~ss_up[.,y_ind];
   
    y_ind = y_ind+1;
	
	if (y_ind == 5)or (y_ind == 7);
		      y_ind = y_ind+1;
	endif;
	
endo;

if mhrun == "3";
   output file= c:\dge_sw\results\shock_sample_A.txt reset; 
elseif mhrun == "4";
    output file= c:\dge_sw\results\shock_sample_B.txt reset; 
elseif mhrun == "5";
     output file= c:\dge_sw\results\shock_sample_C.txt reset; 
endif;

output on;
  
  i = 1;
  do until i > 7;
   
   i "-th   shock  " ;

     "mean    "  "low band   "  "high band"; 
   yy_save[.,1+3*(i-1):3*i];
   i = i + 1;   
   
  endo;
output off; 



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
