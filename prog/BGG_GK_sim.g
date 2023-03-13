/*************************************************************************
**  Simulate the DSGE model 
**************************************************************************/

/* filename: BGG_GK_sim.g
** created:  12/06/2010
** Modifications:
**
*/
new;
library pgraph, user;
  _fcmptol = 1E-10;
  dlibrary -a PszTgSen,PsdTgSen;

#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g

//#include c:\dge_BGG_GK\BGG_GK_mod_1.src
#include c:\dge_BGG_GK\BGG_GK_mod.src

cls;
   
/* Load model parameters
*/


load path = ^priopath lppara[npara,1] = "msim_par.out";
para      = lppara[.,1];   
   /*
   beta    =  0.995 ;
   delta   =  0.025;
   alpha   =  0.3;
   gam_E_ss = 0.972;
   gam_F_ss = 0.972;
   lam      = 0.383;
   psi_w    = 0.05;
   eps      = 11;
   xi_E     = 0.003;
   xi_F     = 0.003;  

   mc_ss    = (eps -1)/eps;
   S_ss     = 1.0075;
   rr_E_ss  = 1.0152;
   rr_F_ss  = rr_E_ss/S_ss;
   rr_ss    = 1/beta;
   nu_ss    = (1-gam_F_ss)*beta*(rr_F_ss-rr_ss)/(1/beta-gam_F_ss); 
   eta_ss   = (1-gam_F_ss)/(1-beta*gam_F_ss);   
   lev_F_ss = eta_ss/(lam-nu_ss);

   K_ss_N_E_ss = 1.919;
   K_ss_Y_ss   = alpha*mc_ss/(rr_E_ss - (1-delta) ); 
   I_ss_Y_ss  = delta* K_ss_Y_ss ; 
   G_ss_Y_ss  = 0.2;
   C_ss_Y_ss  = 1 - I_ss_Y_ss - G_ss_Y_ss ;
   */

   kappa   = para[1];
   h       = para[2];
   sigma_C = para[3];
   sigma_L = para[4];
   phi     = para[5];
   iota_p  = para[6];
   iota_w  = para[7];
   theta_p = para[8];
   theta_w = para[9];
   rho_R   = para[10];
   mu_pi   = para[11];
   mu_y    = para[12];
   rho_a  = para[13];
   rho_c  = para[14];
   rho_k  = para[15];
   rho_e  = para[16];
   rho_f  = para[17]; 
   rho_g  = para[18];
   rho_l  = para[19];  
   e_a    = para[20];
   e_c    = para[21];
   e_e    = para[22];
   e_f    = para[23];
   e_g    = para[24];
   e_k    = para[25];
   e_l    = para[26];
   e_m    = para[27];
   

/* 
** s0:        state of economy in t0
** x0:        initial value of x
** probs0:    agents beliefs about state of economy in t0
** sT:        sequence of states
** yT:        sequence of observed variables
** xT:        sequence of latent variables:
**            [y,inf,r2,xi_y,xi_infl,g,z]'
*/

/* Solve the DSGE model */
/*
** retcode = -1 : non existence                
**         = 0  : existence and uniqueness     
**         = 1  : existence and non-uniqueness
*/

{T1,TC,T0,RC} = dsgesolv(para);    @ rc; T1; TC; T0;  RC;@

@retcode = 0;@

/*
if (RC[1] == 1) AND (RC[2]==1);
   /* determinacy */
   retcode[1] = 0;
   TT = T1;
   RR = TEPS;
   
elseif (RC[1] == 1) AND (RC[2]==0) ;
   /* indeterminacy */
   retcode[1] = 1;
   TT = T1;
   RR = TEPS;  /*
   loglh = loglhzero;
   goto eofevalsoe;
   */
else;
   /* no equilibrium exists, numerical problems */
   retcode[1] = RC[1];
   loglh = loglhzero;
   goto eofsunsim;

endif;
*/     

/* Simulate DSGE model
** Composition of xt: [x(t),pi(t),Et(x(t+1)),Et(pi(t+1)),m(t),deltaM(t),z(t)]'
*/
ZZ      = zeros(11,cols(t1));
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

DD      = zeros(11,1);                              @@

/* Initializations 
*/

Tsimtotal = 200;
Tsimburn  = 200;

yT = zeros(Tsimtotal+Tsimburn,11);                @@
xt = zeros(cols(T1),1);


t = 1;
do until t > Tsimburn + Tsimtotal;

   /* calculate realization of epst */
   epst   = rndn(1,1)*e_a |
            rndn(1,1)*e_c |
            rndn(1,1)*e_e |
            rndn(1,1)*e_f |
            rndn(1,1)*e_g |
            rndn(1,1)*e_k |
            rndn(1,1)*e_l |
            rndn(1,1)*e_m ; 
  
   /* Update vector x(t) */
   xt = T1*xt + TC + T0*epst;   @xt;@
   yT[t,.] = DD' + (ZZ*xt)';

   t = t+1;

endo;

yT = yT[Tsimburn+1:Tsimtotal+Tsimburn,.];

"Data Simulated Successfully";
"Parameter Choice:    " datasel;
"Model Specification: " mspec;
"Parameters (theta):  "; para;
"Mean of YT:          " meanc(yT)';
"StD of YT:           " stdc(yT)';

"RC[1]=" rc[1]  ",  RC[2]=" rc[2]  "     1: success  0: fail ";
/*
"T1=" T1;  
"T0=" T0; 
"TC=" TC;
*/
eofsunsim:


xy(rows(yt), yt);

/********************************************************
** Save simulated data
*/

osimdata         = datapath $+ "\\" $+ lmodel $+ "sim" $+ dataselstr;


create fhsimdata = ^osimdata with SERIES, 11, 8;     @@
wr               = writer(fhsimdata,yT);
closeall fhsimdata;   

"Data saved";
