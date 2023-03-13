/*************************************************************************
**  Simulate the DSGE model 
**************************************************************************/

/* filename: lpsim.g
** created:  10/18/01
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

nirf = 20;
ti    = seqa(1,1,nirf);
   
/* Load model parameters
*/
load path = ^priopath lppara[npara,1] = "msim_par.out";
para      = lppara[.,1];   
   
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
  



{T1,TC,T0,RC} = dsgesolv(para);    @ rc; T1; TC; T0;  RC;@

@retcode = 0;@



/* Simulate DSGE model
** Composition of xt: [x(t),pi(t),Et(x(t+1)),Et(pi(t+1)),m(t),deltaM(t),z(t)]'
*/
ZZ      = zeros(11, rows(T1));                          @@
ZZ[1,1] = 1;             /*  consumption */  
ZZ[2,2] = 1;             /* wage */  
ZZ[3,3] = 1;             /*  labor */  
ZZ[4,4] = 1;              /*  nominal rate  */           
ZZ[5,5] = 1;             /* inflation */  
ZZ[6,7] = 1;              /* output gap */
ZZ[7,9] = 1;             /* External Finance Premium */ 
ZZ[8,18] = 1;            /* Levarage Ratio of Bank */
ZZ[9,19] = 1/5;             /* investment */          
ZZ[10,22] = 1;             /* Levarage Ratio of Firm */  
ZZ[11,23] = 1;           /* Nomial corporate borrowing Rate */     

DD      = zeros(11,1);                              @@

/* Initializations 
*/

yT = zeros(200,11);                @@
xt = zeros(cols(T1),1);


t = 1;
do until t > 200;

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

yT = yT[1:200,.];

"Data Simulated Successfully";
"Parameter Choice:    " datasel;
"Model Specification: " mspec;
"Parameters (theta):  "; para;
"Mean of YT:          " meanc(yT)';
"StD of YT:           " stdc(yT)';

"RC=" rc'  "    1: success  0: fail ";
/*
"T1=" T1;  
"T0=" T0; 
"TC=" TC;
*/
eofsunsim:

sig_chol = eye(8);        @@
//sig_chol[8,8] = 1/4;             @ monetary    @

yyirfall = zeros(nirf,nvar*nshock); 

sh_ind = 1;
do until sh_ind > nshock;

   impact = sig_chol[.,sh_ind];
   yyirf  = zeros(nirf,nvar); nvar;

 
   ss     = (T0*impact);     //rows(zz*ss);
   yyirf[1,.] = (ZZ*ss)';
      
   t = 2;
   do until t > nirf;
      ss = T1*ss;     
      yyirf[t,.] = (ZZ*ss)';
      t = t+1;
   endo;

     yyirfall[.,1+nvar*(sh_ind-1):nvar*sh_ind] = yyirf;

   sh_ind = sh_ind+1;
   
endo;


/*******************************************/

titlestr = "Productivity Shock " $|
           "Preference Shock " $|
           "Firm Net Worth Shock" $|
           "Bank Net wortk Shock" $| 
           "Government Spending Shock " $|
           "Investment Shock " $|
           "Labor Supply Shock " $|          
           "Monetary Policy Shock " ;    @@ 

ystr     =  "Consumption" $| 
            "Real Wage" $|
            "Labor" $|
            "Nominal Rate" $|
             "Inflation" $|
            "Output" $|
            "External Finance Premium" $| 
            "Levarage Ratio of Bank " $| 
              "Investment" $|
            "Levarage Ratio of Firm " $|
            "Nomial corporate borrowing Rate" ;
                       
               


   sh_ind = 1;

      do until sh_ind > 8;
       /* iterate over shocks */

      graphset;
      begwind;
      margin(0,0,0.2,0.2);
      window(3,4,0);           @@
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

       shmult = 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 ;

       y_ind = 1;
       do until y_ind > nvar;
       /* iterate over endogenous variables */
          
          title(titlestr[sh_ind]);
          ylabel(ystr[y_ind]);       

          xy(ti, shmult[sh_ind]*yyirfall[.,nvar*(sh_ind-1)+y_ind]  );

          nextwind;
                 
          y_ind = y_ind+1;
       endo;
       
       sh_ind = sh_ind+1;
    
       endwind; 

     "press any key";    wait;

   endo;