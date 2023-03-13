/* filename:    CEEpm.g
** description: The program maximizes the posterior density of the 
**              DSGE model 
** created:     10/05/2010
*/

new;
 library optmum, pgraph, user; 

 
 dlibrary -a PszTgSen,PsdTgSen;
cls;

format /mat /on /mb1 /ros 16,8;


/* Global defaults
*/
  __output = 1;
  _fcmptol = 1E-10;

#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g
#include c:\dge_BGG_GK\bgg_gk_mod.src;
#include c:\dge_BGG_GK\bgg_gk_lh.src;

/* 
optimization method
*/

opt = 1;   // 1:csminwel  2:optmum

/******************************************************** 
** Import data
** series (nobs,7)
*/

cls;

#include c:\dge_BGG_GK\loaddata.g

nobs  = rows(series_yT);  /* number of observations */
YY    = series_yT;        

   //yy[.,.];  wait;
   
   yy_m = meanc(yy);  
   yy = yy - yy_m';    /* demean */    

    

/******************************************************** 
** Import the information for the prior density
*/
priorfile = lmodel $+ lprior $+ dataselstr $+ subTstr $+ "prt.out";
load path = ^priopath prior[npara,3] = ^priorfile;

pmean  = prior[.,1];
pstdd  = prior[.,2];
pshape = prior[.,3];

/*
pshape = pshape.*para_maskinv[.,mspec];
*/
/******************************************************** 
** Starting Values for Maximization
*/

kappa	=	0.4	;
h	    =	0.6	;
sigma_C	=	1.5;
sigma_L	=	1.3	;
fi   	=	0.04	;
iota_P	=	0.8	;
iota_W	=	0.7	;
theta_P	=	0.8	;
theta_W	=	0.5	;
rho_R	=	0.7	;
mu_Pi	=	1.8	;
mu_Y	=	0.01	;
rho_A	=	0.9	;
rho_C	=	0.8	;
rho_Ak	=	0.9	;
rho_E	=	0.9	;
rho_F	=	0.2	;
rho_G	=	0.9	;
rho_L	=	0.9	;
eA	=	0.4	;
eC	=	1.5	;
eE	=	0.2	;
eF	=	1.0	;
eG	=	1.7;
eK	=	0.4	;
eL	=	2.1	;
eR	=	0.1	;
//uR	=	0.1	;20110407
uC  =   0.3;
uW  =   0.9;
uL  =   0.4;
uPi =   0.9;
uY  =   0.3;
uS  =   0.12;
uLevF=  6.2;
uK  =  0.5;
uLevE=  2.6;
uRE =  0.15;

/*uY	=	0.3	;
uC	=	0.3	;
uK	=	0.5;
uPi	=	0.9	;
uW	=	0.9	;
uL	=	0.4	;
uRE	=	0.15	;
uLevE	=	2.6	;
uLevF	=	6.2;
uS	=	0.12	;*/



para    =  kappa| h | sigma_C | sigma_L | fi |  iota_P |
            iota_W | theta_P | theta_W | rho_R | mu_Pi | mu_Y |
			rho_A | rho_C | rho_AK | rho_E | rho_F | rho_G | rho_L |
            eA | eC | eE | eF | eG | eK |  eL | eR |            
			 uC | uW | uL | uPi | uY | uS | uLevF | uK | uLevE | uRE ;//20110407 uR (down11)|

/*para    =  kappa| h | sigma_C | sigma_L | fi |  iota_P |
            iota_W | theta_P | theta_W | rho_R | mu_Pi | mu_Y |
			rho_A | rho_C | rho_AK | rho_E | rho_F | rho_G | rho_L |
            eA | eC | eE | eF | eG | eK |  eL | eR |            
			 uY | uC | uK | uPi | uW | uL | uRE | uLevE | uLevF | uS ;//20110407 uR (down11)|*/
npara = rows(para);

para_names={kappa h sigma_c sigma_L phi iota_p iota_w  
             theta_p theta_w rho_r  mu_pi  mu_y  
             rho_a   rho_c   rho_k  rho_e  rho_f  rho_g  rho_l
             e_a   e_c  e_e  e_f  e_g e_k e_l e_r
             uC  uW  uL  uPi  uY  uS  uLevF  uK  uLevE  uRE };//20110407 uR (down11) 
/*para_names={kappa h sigma_c sigma_L phi iota_p iota_w  
             theta_p theta_w rho_r  mu_pi  mu_y  
             rho_a   rho_c   rho_k  rho_e  rho_f  rho_g  rho_l
             e_a   e_c  e_e  e_f  e_g e_k e_l e_r
             uY  uC  uK  uPi  uW  uL  uRE  uLevE  uLevF  uS };//20110407 uR (down11)  */

/* Load transformation scheme for parameters
*/
transfile = "trans.out";
load path=^priopath _trspec[npara,4] = ^transfile;


/******************************************************** 
**      Calculate posterior density at starting values
**
*/

cls;
"Prior Density at Starting Value:";  priodens(para,pmean,pstdd,pshape); 

{ lnpY,retcode,obsmean,obsvar } = evaldsge(para,mspec,YY);

/* lnpy = -objfcn(invtrans(para)); */
"Posterior at Starting Value:    " lnpy;
"Press key to continue";
w = keyw;
cls;

/********************************************************
**       Maximize the posterior density
**
*/
cc = 0;
x0 = invtrans(para) ;  


H0= eye(npara)*1E-1;
nit = 1000;
crit= 1E-7;  //-5

/* Use personal optimization procedure
*/

if  opt == 1; 

 {fh,xh,g,H,itct,fcount,retcode} = csminwel(x0,H0,crit,nit,&objfcn);
 
elseif opt == 2;
 
  _opalgr  = 2;       // algorithm; 2:BFGS,  
  _opmiter = 1000;     // maximum number of iteration 
 
  {xh, fout, g, cout}=optmum(&objfcn, x0);

endif;
 

/********************************************************
** Save parameter estimates
*/

paraest = trans(xh);

oparaest         = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pm";
create fhparaest = ^oparaest with PARAEST, 1, 8;
wr               = writer(fhparaest,paraest);
closeall fhparaest;   


/********************************************************
**       Estimation Output
**
*/
{lnpy, retcode, obsmean, obsvar } = evaldsge( paraest,mspec,YY);

lnprio  = priodens( paraest, pmean, pstdd, pshape);

/* Compute gradient in terms of model parameters
*/
gradmod = g./ diag(gradp(&invtrans,paraest));

cls;

output file= c:\dge_BGG_gk\results\post_mode_para.txt reset;

output on;

"Model              " mspec;
"Dataset            " datasel;
"Subsample          " SubT;
"===================================";
"Posterior Mode     " (lnpy+lnprio);
"Likelihood at Mode " lnpy;
"Prior at Mode      " lnprio;
" ";
"Parameter | Estimate |  Start-V  | Prior Mean | Prior Stdd | Gradient";
outmat = para_names'~paraest~trans(x0)~pmean~pstdd~gradmod;
let mask[1,6] = 0 1 1 1 1 1;
let fmt[6,3] =
   "-*.*s " 8 10
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4;

d = printfm(outmat,(0~1~1~1~1~1),fmt);

if lnpy == -1E6 ;
   end;
endif;

output off;

/********************************************************
**       Compute Residuals
*/

/* Moments */

resid1 = YY - obsmean;
resid1 = resid1[2:rows(resid1),.];

resid1_m1 = meanc(resid1)';
resid1_m2 = meanc( (resid1 - resid1_m1)^2 )';
resid1_m3 = meanc( (resid1 - resid1_m1)^3 )';
resid1_m4 = meanc( (resid1 - resid1_m1)^4 )';

resid2 = (YY - obsmean)./sqrt(obsvar);
resid2 = resid2[2:rows(resid2),.];

resid2_m1 = meanc(resid2)';
resid2_m2 = meanc( (resid2 - resid2_m1)^2 )';
resid2_m3 = meanc( (resid2 - resid2_m1)^3 )';
resid2_m4 = meanc( (resid2 - resid2_m1)^4 )';

ti = ti[2:rows(ti)];

nresid = rows(resid1);

/* Autocorrelation */

resid1_ar = zeros(1,7);
resid2_ar = zeros(1,7);

ind = 1;
do until ind > 7;
   resid1_ar[1,ind]=inv(resid1[1:nresid-1,ind]'*resid1[1:nresid-1,ind])
                    *resid1[1:nresid-1,ind]'*resid1[2:nresid,ind];
   resid2_ar[1,ind]=inv(resid2[1:nresid-1,ind]'*resid2[1:nresid-1,ind])
                    *resid2[1:nresid-1,ind]'*resid2[2:nresid,ind];
   ind = ind+1;
endo;

output on;

"";
"";
"Raw Residuals";
"-------------";
"Mean:" resid1_m1;
"Stdd:" sqrt(resid1_m2);
"Skew:" resid1_m3./resid1_m2^(3/2);
"Kurt:" resid1_m4./(3*resid1_m2^2);
"Corr:" resid1_ar;
"";
"Standardized Residuals";
"----------------------";
"Mean:" resid2_m1;
"Stdd:" sqrt(resid2_m2);
"Skew:" resid2_m3./resid2_m2^(3/2);
"Kurt:" resid2_m4./(3*resid2_m2^2);
"Corr:" resid2_ar;
"";

output off;

/****************************
 Drawing Graph of Data  
***************************/
titlestr = "Nominal Interest Rate " $|
           "Output  " $|
           "Consumption " $|
           "Investment " $|
           "Inflation " $|         
           "Wage " $|          
           "Labor " $|
           "Borrowing Rate" $|
           "Bank Leverage " $|
           "Firm Leverage " $|         
           "Spread";    

graphset;
    begwind;
    margin(0,0,0.2,0.2);
    window(4,3,0);           @@
    fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    @_pcolor = {4 4 4};@
    _pltype = {6 1 1};
    _protate = 0; 
    _plctrl = 0;
    _plwidth = 5.5;
    _psymsiz = 5;
    _plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";
       
    y_ind = 1;
    do until y_ind > 11;
       /* iterate over endogenous variables */
          
        title(titlestr[y_ind]);

        if (datasel ==2) or (datasel ==3);
            xtics(1985, 2010,2,4);
            x=seqa(0/4,1/4, rows(yy))+1985;
        endif;
          
        xy(x, yy[.,y_ind] );
        nextwind;
        y_ind = y_ind+1;
    endo;
    endwind;   
end;

/****************************************************/
/*                 PROCEDURES                       */
/****************************************************/

proc (1) = objfcn(para);
local lnpY, lnprio, obsmean, obsvar, shock, modelpara, retcode;

/* likelihood 
*/


modelpara = trans(para);



{lnpY,retcode,obsmean,obsvar} = evaldsge(modelpara,mspec,YY); 



/* Evaluate the Prior distribution
*/

lnprio = priodens(modelpara, pmean, pstdd, pshape);


/*
locate 25,1;
"Likelihood" lnpY  "Return Code" retcode;
*/

retp(real(-lnpY-lnprio));  /* We minize the inverse of the likelihood fcn */
endp;


