
/* filename:    CEEmh_data_rich_C.g
**
** description: The program generates draws from the posterior
**              distribution via the random walk Metropolis Hastings Algorithm
**              
** created:     April/11/2011
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

#include c:\DSGE_NIMN_V0\pathspec.g
#include c:\DSGE_NIMN_V0\BGG_GK_spec.g 

npara = 27;  // number of parameterts

/******************************************************** 
**  Configure the Metropolis Algorithm
*/

file_update = 0;   // 0: start from first draw  1: update from last draw

case = 5;        // 3:data_rich_A,   5:data_rich_B,   0: prior,    2: posterior 
if case == 5;            
   mhrun   = "5_SV";       
elseif  case == 3;   
   mhrun   = "3_SV";
endif;

kk      = 5;         // thining rate 
cc      = 0.2;      // ajustment rate
nblocks = 300;
nsim    = 100;

burnin_blocks = floor(nblocks*0.5);    //  block number of burnin of sampling of state variables

opath = parapath $+ "\\mhrun" $+ mhrun;

_bounds =
    (1E-10 ~ 5      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 5      )|
    (1E-10 ~ 5      )|
    (1E-10 ~ 5      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1     )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 1      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      )|
    (1E-10 ~ 10      );     @@

/******************************************************** 
** Import data
** series (nobs,11)
*/

#include c:\DSGE_NIMN_V0\loaddata.g;   //case_c_zz; wait;   

nobs  = rows(series_yT);  /* number of periods of observations */
YY    = series_yT;   
     
     yy_m = meanc(yy);  
     yy= yy - yy_m';      /* demean */    

ndata  = cols(yy);       /* number of observations */
nstate = 42;             /* number of state variables */

// ZZ = ( ndata times nstate ) matrix 
zz = make_zz(nstate);  //zz; wait;

gam = zz;

ZZ_unobs = zeros(4, nstate); // selection matrix for unobservable states
ZZ_unobs[1,41] = 1;  // Gamma^E: Borrowing rate for Corp
ZZ_unobs[2,42] = 1;  // Gamma^F: Borrowing rate for Bank
ZZ_unobs[3,36] = 1;  // G: Government Expenditure 
ZZ_unobs[4,20] = 4;  // SpreadF


/******************************************************** 
** Import the information for the prior density
*/
priorfile = lmodel $+ lprior $+ dataselstr $+ subTstr $+ "prt.out";
load path = ^priopath prior[npara,3] = ^priorfile;

pmean  = prior[.,1];
pstdd  = prior[.,2];
pshape = prior[.,3];

/********************************* 
** load Posterior Mode Estimates
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
 
// HH = 0.1*eye(ndata);   psi = zeros(ndata,1);
HH = 0.1*eye(ndata);   psi = 0.5*ones(ndata,1);
// posi_def = 0; check the covariance matrix is positive definite

{lnpost0, lnpy0 } = fcn(para,HH, psi,gam);
// {lnpost0, lnpy0 } = fcn(para,HH, psi,gam ,posi_def ); check the covariance matrix is positive definite

"";
"Prior*Likelihood at Mode:" lnpost0;
"Likelihood at Mode      :" lnpy0;
"";
"Press Key to Continue";
k = keyw;
cls;

/******************************************************** 
** load multiplication matrix for parameters
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
//sigpropinv;
//k=keyw;
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
osv   = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "sv";
ogam  = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "gam";

ovariunobs = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "vaunobs";  // output files for unobservable smoothed states

statname = { POST LIKE PROP PROPPOST REJECT };

if file_update == 1;  // update --> yes 
   open fhstat = ^ostat for update;	
   open fhvari = ^ovari  for update;
   open fhpara = ^opara for update;	
   open fhgam  = ^ogam  for update;
   open fhshock = ^oshock for update;
   open fhsv = ^osv for update;	   
   open fhvariunobs =^ovariunobs for update;

	   eofloop = 0;
	   ndraws = 0;
	   do until eofloop; 
   
            drawblock = readr( fhpara, 1000 );
	   
	        ndraws = ndraws + rows(drawblock);
	        eofloop = eof(fhpara);

       endo;
   
else;
     create fhstat=^ostat with ^statname, 5, 8;
     create fhpara=^opara with ^para_names, npara+ndata*3+nshock*4, 8;
     create fhvari=^ovari with ^para_names, (rows(yy)-1)*ndata, 8;
     create fhshock=^oshock with ^para_names, (rows(yy)-1)*nshock, 8;  // number of shock is 8
	 create fhsv=^osv with ^para_names, (rows(yy)-1)*nshock, 8;  // number of shock is 8
     create fhgam  =^ogam with ^para_names, cols(yy)*nstate, 8;
	 
	 create fhvariunobs =^ovariunobs with ^para_names, (rows(yy)-1)*4, 8;

endif;

/******************************************************** 
**              Metropolis Hastings Algorithm
*/

/*  Prior and Save area */
 
 shocksum = zeros( nobs-1, nshock);
 decomp_sum = zeros( nobs-1, nshock*ndata);
 decomp_sum_unobs = zeros( nobs-1, nshock*4 ); // using for hist. decomp. of unobs. three states

/* Initialize PARA NEW
*/

// HH = 0.1*eye(ndata);  psi = 0.1*ones(ndata,1);
HH = 0.1*eye(ndata);  psi = 0.5*ones(ndata,1);

sv_state = 4*ones((nobs-1),nshock,);
sv_para  = zeros(4,nshock); 
           sv_para[1,.]=0.5*ones(1,nshock); sv_para[3,.]=0.5*ones(1,nshock); 
	       sv_para[4,.]=0.5*ones(1,nshock);


valid = 0;
do until valid == 1;

     retcode_1 = 1;
     do until retcode_1 == 0; 
       
          paranew          = para + cc*(sigscale*rndn(npara,1));

          {postnew, likenew} = fcn(paranew, HH, psi, gam);
         //{postnew, likenew} = fcn(paranew, HH, psi, gam, posi_def); check the covariance matrix is positive definite
         
         {T1, TC, TEPS,RC} = dsgesolv(paranew);  /* solve dsge */

         if (RC[1] == 1) and (RC[2]==1);
               /* determinacy */
              retcode_1 = 0;
          elseif (RC[1] == 1) and (RC[2]==0) ;
              /* indeterminacy */
              retcode_1 = 1;
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
 
/* update proc */  
  if file_update == 1;  // update --> yes 
	  drawpara = seekr(fhpara, ndraws); 
	  drawpara = readr(fhpara, 1); 
	  
	  paraold = drawpara[1:npara]'; //paraold'; 
	  sigma = drawpara[npara+1:ndata+npara];    //sigma';
	  HH = zeros(ndata,ndata);
	  
	  iii = 1;
	  do until iii > ndata;
		 HH[iii,iii]=sigma[iii];  
	     iii = iii+1;
	  endo;	 
	  
	  psi = drawpara[npara+ndata+1:2*ndata+npara]'; //psi';
	   
	  drawgam = seekr(fhgam, ndraws); 
	  drawgam = readr(fhgam, 1);  
	  gam   =  reshape(drawgam,nstate,ndata)';  //gam;
	  
	  drawstat = seekr(fhstat, ndraws);
	  drawvari = seekr(fhvari, ndraws);
	  drawshock = seekr(fhshock, ndraws);   
	  drawvariunobs = seekr(fhvariunobs, ndraws);

endif;

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

   sigmasim = zeros( floor(nsim/kk), ndata);
   psisim = zeros( floor(nsim/kk), ndata);
   gamsim = zeros( floor(nsim/kk), ndata);
   gam_mat_sim = zeros(floor(nsim/kk), ndata*nstate);
  
   state_save  = zeros(  floor(nsim/kk), ndata*(nobs-1)  );
   shock_save  = zeros(  floor(nsim/kk), 8*(nobs-1)  );

   state_save_unobs = zeros(floor(nsim/kk), 4*(nobs-1)); // Gamma^E, Gamma^F, G, SpreadF
   
   // stochastic volatility save
   sv_save  = zeros(  floor(nsim/kk), nshock*(nobs-1)  );
   svsim = zeros( floor(nsim/kk), nshock*4);
    
   
   statenew=0;
   shocknew=0;
   decomp_new=0;
   decomp_new_unobs=0;

   if (iblock == 1) and (file_update == 0); //  update --> No 
      i = 1;               
      j = 1;
      postold      = postnew;
      likeold      = likenew;
      paraold      = paranew;     
   else;
      j=1; i=1;
   endif;      
   
   do until j > nsim;
         
    /* Gibbs sampling and smoothing */
	if (j > 1); 
		statenew0 = statenew;
		shocknew0 = shocknew;
		decomp_new0 = decomp_new;
		decomp_new_unobs0 = decomp_new_unobs;
	endif;
	
	{statenew, shocknew, decomp_new, decomp_new_unobs}  = simulation_smoother(gam, HH, paraold, psi, ZZ_unobs);  
	
 	 // stochastic volatility
	{sv_vola, sv_para, sv_state }=  sv_mh(shocknew, sv_state, sv_para ); // SV sample

    {HH}  = sig(statenew, psi, gam);  
    
    {psi } = coef(statenew,HH,gam);  
    
    {gam, gam_g} = coef_g(statenew, HH, psi, gam) ;  
    

    /* MH of Deep parameters */

     paranew          = paraold + cc*(sigscale*rndn(npara,1));

     {postold, likeold } = fcn(paraold, HH, psi, gam);
     {postnew, likenew } = fcn(paranew, HH, psi, gam);



     propdens   = -0.5*SIGPROPDIM*ln(2*pi) - 0.5*SIGPROPlndet - 0.5*SIGPROPDIM*ln(cc^2)
                  -0.5*(paranew - paraold)'*SIGPROPINV*(paranew - paraold)/cc^2;

     r = minc(1 | exp( postnew - postold));  
  
     if rndu(1,1) < r;
       /* Accept proposed jump
       */ 
          paraold = paranew;
          postold = postnew;
          likeold = likenew;

     else;
       /* Reject proposed jump
       */   
         rej[j] = 1;
     endif;        
      
    
   /*  Save  parameter  */ 

     if fmod(j,kk)== 0;         /* thining kk-th sampling */
            propsim[i,1]     = propdens;   
            proppostsim[i,1] = postnew;
            
            psisim[i,.]      = psi';
            gamsim[i,.]      = gam_g';
            gam_mat_sim[i,.] = vec(gam)';
			sigmasim[i,.]    = sqrt(diag(HH)'); 
			svsim[i,.]       = vec(sv_para)';

            likesim[i,1] = likeold;
            postsim[i,1] = postold;
            parasim[i,.] = paraold';
                  
            rej_kk[i] = meanc(rej[j-kk+1:j]); 

            ss = statenew * gam' ;  //cols(vec(ss)'); cols( state_save[i,.]);
			ss_unobs = statenew*ZZ_unobs'; // selecting unobservable states

            state_save[i,.] = vec(ss)'; 
			state_save_unobs[i,.] = (SS_unobs[.,1]|SS_unobs[.,2]|SS_unobs[.,3]|SS_unobs[.,4])'; // unobservable smoothed state (four states)
   
			
            shock_save[i,.] = (shocknew[1,.]~shocknew[2,.]~shocknew[3,.]~shocknew[4,.]~
                          shocknew[5,.]~shocknew[6,.]~shocknew[7,.]~shocknew[8,.]) ;  
						  
			sv_save[i,.]= vec(sv_vola')';			  
                   
             locate 15,1;
            "Block" iblock "of" nblocks;
            "Simulation step" j "of" nsim;
            "Proposal Dens  " propsim[i,1];
            "Parameters     "; parasim[i,.];
            "sigma";           sigmasim[i,.];
            "psi";             psisim[i,.]; 
            "gam";             gamsim[i,.]; 
           
             i = i + 1;  
                        
           if iblock > burnin_blocks;                    
               shocksum = shocksum + shocknew';     
               decomp_sum = decomp_sum + decomp_new; 
			   decomp_sum_unobs = decomp_sum_unobs + decomp_new_unobs; // historical decomp. for unobs. states
		   endif;
           		   
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
    wr = writer( fhpara,parasim~svsim~sigmasim~psisim~gamsim );      
    wr = writer( fhvari, state_save );
    wr = writer( fhshock, shock_save );
    wr = writer( fhgam, gam_mat_sim );
    wr = writer(fhvariunobs, state_save_unobs); // write unobservable states
	wr = writer( fhsv, sv_save );

    tend = date;
    "Elapsed Time"; etstr(ethsec(tstart,tend));
	
	if (iblock > 1+burnin_blocks) and (fmod(iblock,nblocks/10)== 0); 
		#include c:\DSGE_NIMN_V0\hist_decomp_save.src
		#include c:\DSGE_NIMN_V0\hist_decomp_save_unobs.src
//		#include c:\dsge_nimn_v0\plot_shock_unobs.src 
	endif;
  
    
    iblock = iblock+1;
endo;  

 /*  close OUTPUT files */
   closeall fhstat, fhpara, fhvari, fhvariunobs,fhshock, fhgam, fhsv ;  

/* End MH within Gibbs  */

/* plot of shocks and histrical decomp */ 

#include c:\dsge_nimn_v0\plot_shock.src;
#include c:\dsge_nimn_v0\plot_shock_unobs.src 
//#include c:\DSGE_NIMN_V0\hist_decomp_save.src; 


end;


/****************************************************/
/*                 PROCEDURES                       */
/****************************************************/ 

proc (2) = fcn(para, HH, psi, gam);
// proc (2) = fcn(para, HH, psi, gam, posi_def); check the covariance matrix is positive definite
local lnpY, lnprio, obsmean, obsvar, lnpost, retcode, parabd_ind1, parabd_ind2, tt,rr;
 
    parabd_ind1 = para .> _bounds[.,1];
    parabd_ind2 = para .< _bounds[.,2];


    if (parabd_ind1 == 1) and (parabd_ind2 == 1);
//  if (parabd_ind1 == 1) and (parabd_ind2 == 1) and (posi_def == 0); check the covariance matrix is positive definite

       {lnpY,retcode,obsmean,obsvar } = evaldsge_dr(para,mspec,YY, HH, psi,gam ); 
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

#include c:\DSGE_NIMN_V0\bgg_gk_mod.src;
#include c:\DSGE_NIMN_V0\bgg_gk_lh_dr_C.src;  
#include c:\DSGE_NIMN_V0\bgg_gk_data_rich_C.src; 
#include c:\DSGE_NIMN_V0\prog_sv\sv_sampler_leverage.src; 

