/* filename:    varipmom.g
** description: Compute posterior mean and std and CI based on
**              output of posterior simulator
*/

new;
closeall;
library user, pgraph;
cls;
outwidth 128;

#include c:\dge_BGG_GK\prog_case_c\pathspec.g
#include c:\dge_BGG_GK\prog_case_c\BGG_GK_spec.g   


mhrun  = "5";        @ 3:data_rich_A,  4: data_rich_B,  5: data_rich_C   @

nburn       = 10;        /* Number of initial draws to be discarded */
hpdprob     = 0.90;       /*  Credible Band  */  
nblock      = 5;       /*  number of draws for one block */

/*******************************
Load Data
*******************************/

#include c:\dge_BGG_GK\prog_case_c\loaddata.g    

nobs  = rows(series_yT);            /* number of observations */
YY    = series_yT;          //~ series_yt1~ series_yt2~ series_yt3  ;     /* data rich setting */ 
     
     yy_m = meanc(yy);  
     yy= yy - yy_m';      /* demean */ 

ndata = cols(yy);     /* number of data index */     


/******************************************************************
**         Load Parameter Draws from MH Output
*/


lpath  = parapath $+ "\\mhrun" $+ mhrun;
//opath  = parapath $+ "\\mhrun" $+ mhrun;

ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "va";
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
drawrow     = seekr(fhdraw, nburn); 
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
ss_mean = zeros(nobs-1,ndata);
ss_low  = zeros(nobs-1,ndata);
ss_up   = zeros(nobs-1,ndata);

i = 1;   
do until i > ndata;
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

if mhrun == "3" ;
         titlestr = "Consumption " $|
                    "Wage " $|
		            "Labor " $|
                    "Nominal Interest Rate "  $|
            		   "Inflation " $|  
                    "Output  " $|
                    "Spread"  $| 
                    "Bank Leverage " $|
		            "Investment " $|
                    "Firm Leverage " $|                  
                    "Borrowing Rate" ; 
elseif mhrun == "5" ;
         titlestr = "Consumption 1" $|"Consumption 2" $|"Consumption 3" $|"Consumption 4" $|
                    "Wage 1" $| "Wage 2" $| "Wage 3" $| "Wage 4" $| 
		            "Labor 1" $| "Labor 2" $| "Labor 3" $| "Labor 4" $|
                    "Nominal Interest Rate "  $| 
            	    "Inflation 1" $| "Inflation 2" $| "Inflation 3" $| "Inflation 4" $|  
                    "Output 1" $| "Output 2" $| "Output 3" $| "Output 4" $|
                    "Spread 1"  $|  "Spread 2"  $|  "Spread 3"  $|  "Spread 4"  $|
                    "Bank Leverage 1" $| "Bank Leverage 2" $|"Bank Leverage 3" $| "Bank Leverage 4" $|
		            "Investment 1" $| "Investment 2" $| "Investment 3" $| "Investment 4" $|
                    "Firm Leverage 1" $|  "Firm Leverage 2" $| "Firm Leverage 3" $|                   
                    "Borrowing Rate 1" $| "Borrowing Rate 2" $| "Borrowing Rate 3" $| "Borrowing Rate 4"  ; 
endif;

i = 1;
do until i > ndata;
graphset;
    begwind;
    margin(0,0,0.2,0.2);
    window(4,3,0);           @@
    fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    _pcolor = {2 4 3 3};
    _pltype = {6 1 1 1};
    _protate = 0; 
    _plctrl = 0;
    _plwidth = 5.5;
    _psymsiz = 5;
    _plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";
       
       y_ind = 1;

       do until (y_ind > 12) or (i > ndata);
       /* iterate over endogenous variables */
          
          title(titlestr[i]);

          if (datasel ==2 )or (datasel == 3) ;
              xtics(1985, 2010,2,4);
              x=seqa(0/4,1/4, rows(yy)-1)+1985;
          
          endif;
          
          xy(x, yy[2:rows(yy),i]~ss_mean[.,i]~ss_low[.,i] ~ss_up[.,i]);
         
          nextwind;  
                    
          i = i + 1;    
          y_ind = y_ind+1;
       endo;
      
    endwind;   
wait;
 endo;

/* save  data */

 yy_save = {};
 y_ind = 1;

 do until y_ind > ndata;
      yy_save = yy_save~yy[2:rows(yy),y_ind]~ss_mean[.,y_ind]~ss_low[.,y_ind] ~ss_up[.,y_ind];
   
      y_ind = y_ind+1;
    endo;

if mhrun == "3";
       output file= c:\dge_BGG_GK\prog_case_c\results\state_sample_A.txt reset; 

elseif mhrun == "4";
       output file= c:\dge_BGG_GK\prog_case_c\results\state_sample_B.txt reset; 

elseif mhrun == "5";
       output file= c:\dge_BGG_GK\prog_case_c\results\state_sample_C.txt reset; 
endif;

output on;
  
   i = 1;
   do until i > ndata;
   
   i "-th   state variable  " ;

    "       actual_data    "  "mean    "  "low_band   "  "high_band"; 
   yy_save[.,1+4*(i-1):4*i];
   i = i + 1;
   endo;
output off; 

/*
if mhrun == "3";
   output file= c:\dge_BGG_GK\results\state_mean_A.txt reset; 
elseif mhrun == "4";
    output file= c:\dge_BGG_GK\results\state_mean_B.txt reset; 
elseif mhrun == "5";
     output file= c:\dge_BGG_GK\results\state_mean_C.txt reset; 
endif;

output on;
    
    ss_mean;

output off;
*/

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
