/* filename:    parapmom.g
** description: Compute posterior mean and std and CI based on
**              output of posterior simulator
*/

new;
closeall;
library user, pgraph;
cls;
outwidth 128;

mhrun  = "5";        @ 3:case A,   0: prior,    2:regular DGE    @

nburn       = 10;            /* Number of initial draws to be discarded */
hpdprob     = 0.90;         /*  Credible Band  */  
nblock      = 5;        /*  number of draws for one block */

/******************************************************************
**         Load Parameter Draws from MH Output
*/
#include c:\dge_BGG_GK\prog_case_c\pathspec.g
#include c:\dge_BGG_GK\prog_case_c\BGG_GK_spec.g

lpath  = parapath $+ "\\mhrun" $+ mhrun;
opath  = parapath $+ "\\mhrun" $+ mhrun;

ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa";
odraws = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pas";

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

drawmean   = drawmean/ndraws; //check ndraws 2011/04/13
drawsqmean = drawsqmean/ndraws;
drawstdd   = sqrt(drawsqmean - (drawmean)^2);
drawcov    = drawcross/ndraws - drawmean'*drawmean;

/* Report Posterior Mean and stderror
**
*/
"MODEL       " lmodel;
"PRIOR       " lprior;
"MH-Run      " mhrun;
"Subsample   " subT;
"Dataset     " dataselstr;

"=============================================";
" " ;
 

para_names={kappa h sigma_c sigma_L phi iota_p iota_w  
             theta_p theta_w rho_r  mu_pi  mu_y  
             rho_a   rho_c   rho_k  rho_e  rho_f  rho_g  rho_l
             e_a   e_c  e_e  e_f  e_g e_k e_l e_r };


"Parameter | Mean  |  StdD |";

outmat = para_names'~drawmean[.,1:27]' ~ drawstdd[.,1:27]' ; 

let mask[1,3] = 0 1 1;
let fmt[3,3] =
   "-*.*s " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4;

d = printfm(outmat,mask,fmt);


/* Part 2: HPD Interval
*/
@cls;@
drawrow     = seekr( fhdraw, nburn); 
drawci      = zeros(2,drawdim);
drawcd      = zeros(3,drawdim);

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
      locate 1,1;
      "Part 2";
      "ndraws" ndraws;
      "Column" j;
         
   endo;

   
   drawcd[.,j] = conv_diag(drawcol,ndraws);
   
   drawci[.,j] = hpdint(drawcol,hpdprob);
   j = j+1;

endo;
clear drawcol;

/* Report Posterior Mean and stderror
**
*/

outmat = para_names'~real(drawmean[.,1:27]')~ real(drawcd[3,1:27])'
        ~ real(drawstdd[.,1:27]')~ real(drawci[.,1:27]') ~ real(drawcd[1:2,1:27])';

draw_names_A = { StD_1 StD_2 StD_3 StD_4 StD_5 StD_6 StD_7 StD_8 StD_9 StD_10 StD_11
                 rho_1 rho_2 rho_3 rho_4 rho_5 rho_6 rho_7 rho_8 rho_9 rho_10 rho_11  };

draw_names_C = { StD_1 StD_2 StD_3 StD_4 StD_5 StD_6 StD_7 StD_8 StD_9 StD_10
                 StD_11 StD_12 StD_13 StD_14 StD_15 StD_16 StD_17 StD_18 StD_19 StD_20
                 StD_21 StD_22 StD_23 StD_24 StD_25 StD_26 StD_27 StD_28 StD_29 StD_30
                 StD_31 StD_32 StD_33 StD_34 StD_35 StD_36 StD_37 StD_38 StD_39 StD_40
                 rho_1 rho_2 rho_3 rho_4 rho_5 rho_6 rho_7 rho_8 rho_9 rho_10
                 rho_11 rho_12 rho_13 rho_14 rho_15 rho_16 rho_17 rho_18 rho_19 rho_20
                 rho_21 rho_22 rho_23 rho_24 rho_25 rho_26 rho_27 rho_28 rho_29 rho_30
                 rho_31 rho_32 rho_33 rho_34 rho_35 rho_36 rho_37 rho_38 rho_39 rho_40
                 Gam_1  Gam_2 Gam_3 Gam_4 Gam_5 Gam_6 Gam_7 Gam_8 Gam_9 Gam_10
                 Gam_11  Gam_12 Gam_13 Gam_14 Gam_15 Gam_16 Gam_17 Gam_18 Gam_19 Gam_20
                 Gam_21  Gam_22 Gam_23 Gam_24 Gam_25 Gam_26 Gam_27 Gam_28 Gam_29 Gam_30
                 Gam_31  Gam_32 Gam_33 Gam_34 Gam_35 Gam_36 Gam_37 Gam_38 Gam_39 Gam_40  };


if    mhrun  == "3";
    outmat = outmat|
             draw_names_A'~real(drawmean[.,28:49]')~ real(drawcd[3,28:49]')~ real(drawstdd[.,28:49]')
             ~ real(drawci[.,28:49]')~ real(drawcd[1:2,28:49])' ;

elseif  mhrun  == "5";
    outmat = outmat|
             draw_names_C'~real(drawmean[.,28:147]')~ real(drawcd[3,28:147]')~ real(drawstdd[.,28:147]')
             ~ real(drawci[.,28:147]')~ real(drawcd[1:2,28:147])' ;


/*
elseif mhrun == "4" ;
     outmat = outmat|
             draw_names_B'~real(drawmean[.,28:62]')~real(drawcd[3,28:62]')~real(drawstdd[.,28:62]')
             ~ real(drawci[.,28:62]')~ real(drawcd[1:2,28:62]') ;
*/
endif; 




 output on;

cls;
"MODEL       " lmodel;
"PRIOR       " lprior;
"MH-Run      " mhrun;
"Subsample   " subT;
"Coverage    " hpdprob;  
"Dataset     " dataselstr;
"burnin      " nburn;
"Draws       " ndraws;
"=============================================";
" " ;


"Parameter | Mean  |   SE   |  StdD | CI(LOW) | CI(HIGH)|  CD  | inefficiency factor";

let mask[1,8] = 0 1 1 1 1 1 1 1;
let fmt[8,3] =
   "-*.*s " 8 10
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4;

d = printfm(outmat,mask,fmt);


create fhdraws = ^odraws with DRAWSUM, drawdim, 8;
wr = writer(fhdraws, real(drawmean | drawstdd | drawcov | drawci) );
closeall fhdraws;

 output off;

/***********************/
/*
if mhrun == "3";
     output file= c:\dge_bgg_gk\results\para_a.txt reset; 
elseif mhrun == "4";
     output file= c:\dge_bgg_gk\results\para_b.txt reset;
elseif mhrun == "5";
     output file= c:\dge_bgg_gk\results\para_c.txt reset;
endif; 


 output on;

   real(drawmean[1:27]);
  
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

/****************************************************/

proc (1) = conv_diag(drawcol,ndraws);

local   cd, a, b, z, k, n, w, mean_a, mean_b, mean_n, i, 
        sigma_a, sigma_b, sigma_n, 
        sigma_a_sum, sigma_b_sum, sigma_n_sum;

cd = zeros(3,1);

//  Geweke convergence diagonostic 

w = round(ndraws*0.1);
A = round(ndraws*0.1);
B = round(ndraws*0.5);

mean_A = meanc(drawcol[1:A]);
mean_B = meanc(drawcol[B:ndraws]);

sigma_A = zeros(w+1,1);
sigma_B = zeros(w+1,1);

sigma_a_sum = 0; 
sigma_b_sum = 0;

i = 0;
do until i > w;
   sigma_A[i+1] = (drawcol[1:A]-mean_A)'*(drawcol[1+i:A+i]-mean_A)/A^2;
   
   sigma_B[i+1] = (drawcol[B:ndraws]-mean_B)'*(drawcol[B-i:ndraws-i]-mean_B)/B^2;

   // parzen window
   if i > 0;
      z = i / w;
      
      if  z < 1/2;
             k = 1 - 6*z^2 + 6*z^3;
      else;
             k = 2*(1-z)^3;
      endif;  

      sigma_A_sum  = sigma_A_sum + (2 * A)/(A-1)*k*sigma_A[i+1];

      sigma_B_sum  = sigma_B_sum + (2 * B )/(B-1)*k*sigma_B[i+1];

   endif;

   i = i + 1 ; 
endo;

cd[1,1] = (mean_A - mean_B)/sqrt(sigma_A[1]+sigma_B[1]+sigma_a_sum+sigma_b_sum);


//  inefficiency factor

w = round(ndraws*0.1);
n = ndraws;

mean_n = meanc(drawcol[1:n-w]);

sigma_n = zeros(w+1,1);

sigma_n_sum = 0; 

i = 0;
do until i > w;
   sigma_n[i+1] = (drawcol[1:n-w]-mean_n)'*(drawcol[1+i:n-w+i]-mean_n)/n^2;
   
   // parzen window
   if i > 0;
      z = i / w;
      
      if  z < 1/2;
             k = 1 - 6*z^2 + 6*z^3;
      else;
             k = 2*(1-z)^3;
      endif;  

      sigma_n_sum  = sigma_n_sum + (2 * n)/(n-1)*k*sigma_n[i+1];
      
   endif;

   i = i + 1 ; 
endo;

cd[2,1] = 1 + sigma_n_sum/sigma_n[1] ;  // inefficiency factor

cd[3,1] = sqrt( sigma_n[1] + sigma_n_sum );  // standard error

retp(cd);
endp
