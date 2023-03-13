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

#include c:\dge_sw\pathspec.g    
#include c:\dge_sw\ceespec.g    


/******************************************************** 
**  Configure the Metropolis Algorithm
*/
mhrun   = "3";      // 3:data_rich_A,   0: prior,    2: posterior   
kk   =    10;         // thining rate 
cc      = 0.2;      // ajustment rate
nblocks =  300;
burnin_blocks = 1;    //  block number of burnin of sampling of state variables
nsim    = 1000;

opath = parapath $+ "\\mhrun" $+ mhrun;

_bounds =
    (1E-5 ~ 1      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1e-5 ~ 10      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1     )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 5      )|
    (-1  ~ 5      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 1      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 50      )|
    (1E-5 ~ 10    )|
    (1E-5 ~ 50      )|
    (1E-5 ~ 10      )|
    (1E-5 ~ 10      );     @@

/******************************************************** 
** Import data
** series (nobs,7)
*/
#include c:\dge_sw\loaddata.g    

nobs  = rows(series_yT);  /* number of observations */
YY    = series_yT;   
     
     yy_m = meanc(yy);  
     yy= yy - yy_m';      /* demean */    


nstate = 21;               /* number of state variables */
ZZ   = zeros(7,nstate);
ZZ[1,1] = 1;               /* output gap */
ZZ[2,2] = 4;               /* inflation */
ZZ[3,3] = 1;               /* wage */
ZZ[4,6] = 1;               /* investment */            
ZZ[5,7] = 1;               /*  consumption */
ZZ[6,8] = 4;               /*  nominal rate  */
ZZ[7,10] = 1;              /*  labor */  
      
/******************************************************** 
** Import the information for the prior density
*/
priorfile = lmodel $+ lprior $+ dataselstr $+ subTstr $+ "prt_dr.out";
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

print "mode=" para;

"Press Key to Continue";
wait;

/******************************/

//load path = ^priopath para_names;

para_names={h  sigma_c  sigma_L  phi  phi0  psi gam_p 
            gam_w  xi_p   xi_w  rho_m 
            mu_pi  mu_y  rho_z  rho_c  rho_g  rho_L  rho_i 
            e_c  e_inv  e_q e_L   e_w   e_z  e_p  e_g  e_m };

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
 
HH = 0.0*eye(7);   psi = zeros(7,1);
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

lhess = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "hs";
open fhhess = ^lhess for read;
SIGPROPINV  = readr( fhhess,npara );
closeall fhhess;

 
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
create fhpara=^opara with ^para_names, npara+7*2, 8;
create fhvari=^ovari with ^para_names, (rows(yy)-1)*7, 8;
create fhshock=^oshock with ^para_names, (rows(yy)-1)*9, 8;

/******************************************************** 
**              Metropolis Hastings Algorithm
*/


/*  Prior and Save area */

 R0_V=EYE(1)/1;
 T0_V=0;            /*  PSI; autocorrelations of measurement errors */
 
 D0_=0;
 V0_=0;             /*  Sigma^2;  variances of measurement errors */

 
 shocksum = zeros( nobs-1, nshock);
 decomp_sum = zeros( nobs-1, nshock*7);

/* Initialize PARA NEW
*/

HH = 0.01*eye(7);  psi = 0.5*ones(7,1);

valid = 0;
do until valid == 1;

   retcode_1 = 1;
     do until retcode_1 == 0; 
       
          paranew          = para + cc*(sigscale*rndn(npara,1));

          paranew[25]=0.0001; paranew[23]=0.0001; //paranew[21]=0.0001;


         {postnew, likenew} = fcn(paranew, HH, psi );
 
         {T1, TC, TEPS,RC} = dsgesolv(para);  /* solve dsge */


          if (RC[1] == 1) and (RC[2]==1);
             /* determinacy */
             retcode_1 = 0;
          endif;

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

   sigmasim = zeros( floor(nsim/kk), 7);
   psisim = zeros( floor(nsim/kk), 7);
  
   state_save  = zeros(  floor(nsim/kk), 7*(nobs-1)  );
   shock_save  = zeros(  floor(nsim/kk), 9*(nobs-1)  );

   
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

       paranew[25]=0.0001; paranew[23]=0.0001; //paranew[21]=0.0001;


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

            state_save[i,.] = (ss[., 1]|ss[., 2]|ss[., 3]|ss[., 4]|ss[., 5]|ss[., 6]|ss[., 7])';     
            shock_save[i,.] = (shocknew[1,.]~shocknew[2,.]~shocknew[3,.]~shocknew[4,.]~
                               shocknew[5,.]~shocknew[6,.]~shocknew[7,.]~shocknew[8,.]~shocknew[9,.]) ;  
                   
          
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

#include c:\dge_sw\prog_dr\plot_shock.src; 

    
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

#include c:\dge_sw\ceemod.src;   
#include c:\dge_sw\prog_dr\ceelh_dr_A.src;  
#include c:\dge_sw\prog_dr\data_rich_A.src; 

