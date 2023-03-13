
/******************************************************************/
/*                                                                */
/*              Plot the impulse response functions               */
/*                                                                */
/*                        Frank Schorfheide                       */
/*                                                                */
/******************************************************************/

/* filename: irfplot.g
** created:  09/27/01
*/
new;
library pgraph;
dlibrary -a PszTgSen,PsdTgSen;

cls;

#include c:\dge_sw\pathspec.g   
#include c:\dge_sw\ceespec.g    

nirf  = 20;                /* number of impulse responses */   @@
ti    = seqa(1,1,nirf);

/* Read Posterior Simulation Summary Files
*/
mhrun  = "2";    @  2: SW model, 3: Case A, 4: Case B, 5: Case C   @

lpath  = resupath $+ "\\mhrun" $+ mhrun;
ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "irs";

open fhdraw = ^ldraw for read;

irf_m = readr(fhdraw,1);
irf_s = readr(fhdraw,1);
irf_hpd_l = readr(fhdraw,1);
irf_hpd_h = readr(fhdraw,1);

irf_m     = reshape(irf_m,nvar*nshock,nirf)';
irf_hpd_l = reshape(irf_hpd_l,nvar*nshock,nirf)';
irf_hpd_h = reshape(irf_hpd_h,nvar*nshock,nirf)';

titlestr = "Preference Shock " $|
           "Investment Shock " $|
           "Equity Premium Shock " $|
           "Labor Supply Shock " $|
           "Wage Markup Shock " $|   
           "Productivity Shock " $|
           "Price Markup Shock " $|                     
           "Government Spending Shock " $|
           "Monetary Policy Shock " ;    @@ 

ystr     = "Output" $|
            "Consumption" $|
           "Inflation" $|
             "Nominal Rate" $|
           "Real Wage" $|
            "Rental Rate" $|
             "Captal Stock"  $|  
           "Investment" $|             
           "Labor";     @@

shmult = 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 ;

   
/*******************************************/
   sh_ind = 1;

      do until sh_ind > 9;
       /* iterate over shocks */

      graphset;
      begwind;
      margin(0,0,0.2,0.2);
      window(3,3,0);           @@
      fonts("microb");
      _ptitlht = 0.4;
      _paxht = 0.4;
      _pnumht = 0.35;
      _pcolor = {10 10 10};
      _pltype = {6 1 1};
      _protate = 0; 
      _plctrl = 0;
      _plwidth = 5.5;
      _psymsiz = 5;
      _plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";

       
       y_ind = 1;
       do until y_ind > nvar;
       /* iterate over endogenous variables */
          
          title(titlestr[sh_ind]);
          ylabel(ystr[y_ind]);       

          xy(ti, shmult[sh_ind]*irf_m[.,nvar*(sh_ind-1)+y_ind] 
                 ~ shmult[sh_ind]*irf_hpd_l[.,nvar*(sh_ind-1)+y_ind] 
                 ~ shmult[sh_ind]*irf_hpd_h[.,nvar*(sh_ind-1)+y_ind] );

          nextwind;
                 
          y_ind = y_ind+1;
       endo;
       
       sh_ind = sh_ind+1;
    
       endwind; 

     "press any key";    wait;

   endo;
/********************************/


         sh_ind =1;   
      file = zeros(20, 3);

   if mhrun  == "2";
      OUTPUT FILE=c:\dge_sw\results\irf_sw.txt RESET;
   elseif mhrun == "3";
       OUTPUT FILE=c:\dge_sw\results\irf_A.txt RESET;
   elseif mhrun == "4";
       OUTPUT FILE=c:\dge_sw\results\irf_B.txt RESET;
   elseif mhrun == "5";
       OUTPUT FILE=c:\dge_sw\results\irf_C.txt RESET;
   endif;
 
   do until  sh_ind > 9;
   
       y_ind = 1;
      do until  y_ind > nvar;
      file[.,.] = (shmult[sh_ind]*irf_m[.,nvar*(sh_ind-1)+y_ind]*100 
                 ~ shmult[sh_ind]*irf_hpd_l[.,nvar*(sh_ind-1)+y_ind]*100 
                 ~ shmult[sh_ind]*irf_hpd_h[.,nvar*(sh_ind-1)+y_ind]*100) ;
 

      "    ";
      "impulse response of " ystr[y_ind] "  to  " titlestr[sh_ind] ;
      "   |      mean      |       low     |      upper     |  "
      file ;

      y_ind = y_ind + 1;
      endo;

     sh_ind = sh_ind + 1;

  endo;
          
  OUTPUT OFF;