/* filename:    CEEmh_data_rich_A.g
**
** description: The program generates draws from the posterior
**              distribution via the random walk Metropolis Hastings Algorithm
**              
** created:     10/4/2010
** modified:   
*/

new;
closeall;
library user, pgraph;
 dlibrary -a PszTgSen,PsdTgSen;

cls;

/* Global defaults 
*/
  __output = 1;
  _fcmptol = 1E-10;

#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g 

npara = 27;  // number of parameterts

/******************************************************** 
**  Configure the Metropolis Algorithm
*/
mhrun   = "3";      // 3:data_rich_A,   0: prior,    2: posterior   
kk   =    1;         // thining rate 
cc      = 0.2;      // ajustment rate
nblocks =  10;
burnin_blocks = 1;    //  block number of burnin of sampling of state variables
nsim    = 100;

opath = parapath $+ "\\mhrun" $+ mhrun;

_bounds =
    (1E-5 ~ 5      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 5      )|
    (1E-5 ~ 5      )|
    (1E-5 ~ 5      )|
    (1e-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1     )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 30      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      );     @@

/******************************************************** 
** Import data
** series (nobs,11)
*/
#include c:\dge_BGG_GK\loaddata.g   

nobs  = rows(series_yT);  /* number of observations */
YY    = series_yT;   
     
     yy_m = meanc(yy);  
     yy= yy - yy_m';      /* demean */    


nstate = 42;               /* number of state variables */
ZZ      = zeros(11,nstate);
ZZ[1,1] = 1;               /*  consumption */
ZZ[2,2] = 1;               /* wage */
ZZ[3,3] = 1;              /*  labor */  
ZZ[4,4] = 1;               /*  nominal rate  */
ZZ[5,5] = 1;               /* inflation */
ZZ[6,7] = 1;               /* output gap */
ZZ[7,9] = 1;              /*  spread */
ZZ[8,18] = 1;               /*  bank's leverage  */
ZZ[9,19] = 1;               /* investment */            
ZZ[10,22] = 1;             /*  Firms's leverage */
ZZ[11,23] = 1;              /*  Firms borrowing rate */ 
 
 
      
/******************************************************** 
** Import the information for the prior density
*/
priorfile = lmodel $+ lprior $+ dataselstr $+ subTstr $+ "prt.out";
load path = ^priopath prior[npara,3] = ^priorfile;

pmean  = prior[.,1];
pstdd  = prior[.,2];
pshape = prior[.,3];

/********************************* 
** Load Posterior Mode Estimates
*/

lpara = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pm";
open fhpara = ^lpara for read;
para = readr(fhpara,npara);
closeall fhpara;

//para = para[1:27];  //HH = para[28:38];

print "mode=" para;

"Press Key to Continue";
wait;

/******************************/

//load path = ^priopath para_names;

para_names={kappa h sigma_c sigma_L phi iota_p iota_w  
             theta_p theta_w rho_r  mu_pi  mu_y  
             rho_a   rho_c   rho_k  rho_e  rho_f  rho_g  rho_l
             e_a   e_c  e_e  e_f  e_g e_k e_l e_r };
             //uR  uY  uC  uK  uPi  uW  uL  uRE  uLevE  uLevF  uS };

"Model            " $lmodel;
"Data set         " $dataselstr;
"Prior            " $lprior;
"Subsample        " $subTstr;
"MHrun Index      " $mhrun;
"Scaling of Hess  " cc;
"";
"===================================";
"";
"";
"Parameter   | Estimate ";
outmat = para_names'~para;
let mask[1,2] = 0 1;
let fmt[2,3] =
   "-*.*s " 10 10
   "*.*lf " 10 4;
d = printfm(outmat,(0 ~ 1),fmt);
 
HH = 0.1*eye(11);   psi = zeros(11,1);

{lnpost0, lnpy0 } = fcn(para,HH, psi  );

"";
"Prior*Likelihood at Mode:" lnpost0;
"Likelihood at Mode      :" lnpy0;
"";
"Press Key to Continue";
k = keyw;
cls;



/******************************************************** 
** Load multiplication matrix for parameters
*/

lmult = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "mu";
open fhmult = ^lmult for read;
SIGSCALE   = readr( fhmult,npara );
closeall fhmult;

SIGSCALE =  SIGSCALE[1:27,1:27];


lhess = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "hs";
open fhhess = ^lhess for read;
SIGPROPINV  = readr( fhhess,npara );
closeall fhhess;

sigpropinv = sigpropinv[1:27,1:27];  
 
SIGPROPdim = sumc(ones(npara,1));
{u , s, v} = svd1(SIGPROPINV);    
SIGPROPlndet = 0;

i = 1;
do until i > npara;
   if i <= SIGPROPdim;
      s[i,i]    = 1/s[i,i];
      SIGPROPlndet = SIGPROPlndet + ln(s[i,i]);
   endif;
   i = i+1;
endo;

/******************************************************** 
** Initialize Output files
*/
ostat = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "st";
opara = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa";
ovari = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "va";
oshock= opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "sck";

statname = { POST LIKE PROP PROPPOST REJECT };

create fhstat=^ostat with ^statname, 5, 8;
create fhpara=^opara with ^para_names, npara+11*2, 8;
create fhvari=^ovari with ^para_names, (rows(yy)-1)*11, 8;
create fhshock=^oshock with ^para_names, (rows(yy)-1)*8, 8;

/******************************************************** 
**              Metropolis Hastings Algorithm
*/


/*  Prior and Save area */

 R0_V=EYE(1)/1;
 T0_V=0;            /*  PSI; autocorrelations of measurement errors */
 
 D0_=0;
 V0_=0;             /*  Sigma^2;  variances of measurement errors */

 
 shocksum = zeros( nobs-1, nshock);
 decomp_sum = zeros( nobs-1, nshock*11);

/* Initialize PARA NEW
*/

HH = 0.1*eye(11);  psi = 0.1*ones(11,1);

valid = 0;
do until valid == 1;

   retcode_1 = 1;
     do until retcode_1 == 0; 
       
          paranew          = para + cc*(sigscale*rndn(npara,1));

          {postnew, likenew} = fcn(paranew, HH, psi );
 
         {T1, TC, TEPS,RC} = dsgesolv(para);  /* solve dsge */

         if (RC[1] == 1) and (RC[2]==1);
               /* determinacy */
              retcode_1 = 0;
          elseif (RC[1] == 1) and (RC[2]==0) ;
              /* indeterminacy */
              retcode_1 = 0;
          endif;

         // retcode_1 ;  RC;  postnew; likenew;  wait;  

     endo;

   /* proposal density  */
   propdens  = -0.5*SIGPROPDIM*ln(2*pi) - 0.5*SIGPROPlndet - 0.5*SIGPROPDIM*ln(cc^2)
               -0.5*(paranew - para)'*SIGPROPINV*(paranew - para)/cc^2;

   if postnew > -1E+20;
      valid = 1;
   endif;

endo;   
 


/* Start MH with in Gibbs  */

iblock = 1;
do until iblock > nblocks;

   tstart = date;

   parasim = zeros( floor(nsim/kk),npara);
   likesim = zeros( floor(nsim/kk) ,1);
   postsim = zeros( floor(nsim/kk), 1);
   rej     = zeros( nsim ,1);
   rej_kk  = zeros( floor(nsim/kk) ,1);
   propsim = zeros( floor(nsim/kk) ,1);
   proppostsim = zeros( floor(nsim/kk) ,1);

   sigmasim = zeros( floor(nsim/kk), 11);
   psisim = zeros( floor(nsim/kk), 11);
  
   state_save  = zeros(  floor(nsim/kk), 11*(nobs-1)  );
   shock_save  = zeros(  floor(nsim/kk), 8*(nobs-1)  );

   
   if iblock == 1;
      i = 1;               
      j = 2;
      postold      = postnew;
      likeold      = likenew;
      paraold      = paranew;     
   else;
      j=1; i=1;
   endif;      

   do until j > nsim;

         
    /* Gibbs sampling and smoothing */

    {statenew, shocknew, decomp_new}  = simulation_smoother(T1, TEPS, ZZ, HH, paranew, psi );  

    {HH}  = sig(statenew, psi);  

    {psi } = coef(statenew,HH);  


    /* MH of Deep parameters */

     paranew          = paraold + cc*(sigscale*rndn(npara,1));

     {postold, likeold } = fcn(paraold, HH, psi );

     {postnew, likenew } = fcn(paranew, HH, psi );

     propdens   = -0.5*SIGPROPDIM*ln(2*pi) - 0.5*SIGPROPlndet - 0.5*SIGPROPDIM*ln(cc^2)
                  -0.5*(paranew - paraold)'*SIGPROPINV*(paranew - paraold)/cc^2;

    
   /*  Save  parameter  */ 

     if fmod(j,kk)==0;         /* thining kk-th sampling */
            propsim[i,1]     = propdens;
            proppostsim[i,1] = postnew;
            psisim[i,.]      = psi';
            sigmasim[i,.]    = sqrt(diag(HH)');   

            ss = statenew * ZZ' ;

       state_save[i,.] = (ss[., 1]|ss[., 2]|ss[., 3]|ss[., 4]|ss[., 5]|ss[., 6]
                         |ss[., 7]|ss[., 8]|ss[., 9]|ss[., 10]|ss[., 11])';
     
       shock_save[i,.] = (shocknew[1,.]~shocknew[2,.]~shocknew[3,.]~shocknew[4,.]~
                          shocknew[5,.]~shocknew[6,.]~shocknew[7,.]~shocknew[8,.]) ;  
                   
          
     endif;  
         

     r = minc(1 | exp( postnew - postold));  
  
    if rndu(1,1) < r;

       /* Accept proposed jump
       */
       if fmod(j,kk)==0;             /* thining kk-th sampling */
           postsim[i,1] = postnew;
           likesim[i,1] = likenew;
           parasim[i,.] = paranew'; 
                     
       endif;

          paraold = paranew;
          postold = postnew;
          likeold = likenew;

     else;
       /* Reject proposed jump
       */   
         rej[j] = 1;

         if fmod(j,kk)==0;         /* thining kk-th sampling */
            likesim[i,1] = likeold;
            postsim[i,1] = postold;
            parasim[i,.] = paraold';
                  
            rej_kk[i] = 1;
         endif;
       
     endif;
     
     if fmod(j,kk)==0;
            locate 10,1;
            "Block" iblock "of" nblocks;
            "Simulation step" j "of" nsim;
            "Proposal Dens  " propsim[i,1];
            "Parameters     "; parasim[i,.];
            "sigma";           sigmasim[i,.];
            "psi";             psisim[i,.]; 

           if iblock > burnin_blocks;                    
               shocksum = shocksum + shocknew';     
               decomp_sum = decomp_sum + decomp_new; 
           endif;

            i = i + 1;  
     endif;           

     
     j = j+1;
    endo;

   locate 1,1;
   "Block" iblock "of" nblocks;
   "Rejection: perct " meanc(rej);
   "Likelihood       " meanc(likesim);
   "Posterior        " meanc(postsim);
   "Parameters       "; meanc(parasim)'; 

    wr = writer( fhstat,postsim~likesim~propsim~proppostsim~rej_kk );
    wr = writer( fhpara,parasim~sigmasim~psisim );      
    wr = writer( fhvari, state_save );
    wr = writer( fhshock, shock_save );


    tend = date;
    "Elapsed Time"; etstr(ethsec(tstart,tend));
  
    
    iblock = iblock+1;
endo;  

 /*  close OUTPUT files */
   closeall fhstat, fhpara, fhvari, fhshock ;  

/* End MH within Gibbs  */


/* plot of shocks and histrical decomp */ 

#include c:\dge_BGG_GK\prog_dr\plot_shock.src; 

    
end;


/****************************************************/
/*                 PROCEDURES                       */
/****************************************************/ 

proc (2) = fcn(para, HH, psi);
local lnpY, lnprio, obsmean, obsvar, lnpost, retcode, parabd_ind1, parabd_ind2;
 
    parabd_ind1 = para .> _bounds[.,1];
    parabd_ind2 = para .< _bounds[.,2];


    if (parabd_ind1 == 1) and (parabd_ind2 == 1);

       {lnpY,retcode,obsmean,obsvar } = evaldsge_dr(para,mspec,YY, HH, psi ); 
       lnpy    = real(lnpy);
       lnprio  = priodens( para, pmean, pstdd, pshape);
       lnpost  = lnpy + real(lnprio);

    else;
      
       lnpost = -1E20;
       lnpy   = -1E20;

         locate 20,1;
         "Out of Band of parameters" ;

    endif;
    
retp(lnpost,lnpy);  
endp;

/******************************************************/

#include c:\dge_BGG_GK\bgg_gk_mod.src;
#include c:\dge_BGG_GK\prog_dr\bgg_gk_lh_dr_A.src;  
#include c:\dge_BGG_GK\prog_dr\bgg_gk_data_rich_A.src; 

