/* filename:    BGG_GK_pri1.g
** description: The program generates draws from the prior of the 
**              DSGE model
** created:     12/07/2010
*/

library pgraph, user;
                         cls;  
block_n    = 10;
block_size = 500;
mhrun      = "0";   /* NOTE: we use mhrun = 0 to denote prior */


#include c:\dge_BGG_GK\pathspec.g
#include c:\dge_BGG_GK\BGG_GK_spec.g

opath      = parapath $+ "\\mhrun" $+ mhrun;

/******************************************************** 
** Import the information for the prior density
*/
priorfile = lmodel $+ lprior $+ dataselstr $+ subTstr $+ "prt.out";
load path = ^priopath prior[npara,3] = ^priorfile;

pmean  = prior[.,1];
pstdd  = prior[.,2];
pshape = prior[.,3];



/******************************************************** 
**  Initialize Output files
*/

otheta = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pa0";
create fhtheta = ^otheta with THETA, npara, 8;

/* Loop to generate parameter draws
*/
block_it = 1;
para_acpt = 0;

do until block_it > block_n;

   /* Generate THETA draws
   */

   parasim_t = zeros(block_size,npara);

   i = 1; 
   do until i > npara;
   
      if pshape[i] == 1;    /* BETA Prior */

          j = 1;
        do until j > block_size;
             a = (1-pmean[i])*pmean[i]^2/pstdd[i]^2 - pmean[i];
             b = a*(1/pmean[i] - 1);
      
          if i > 10;    flag = 0;
                do until  flag > 1; 
                       parasim_t[j,i] = rndbeta(1,1,a,b);
                    if parasim_t[j,i] < 0.95;
                         flag = 2;  
                   endif;
                endo; 
          else;  
                  parasim_t[j,i] = rndbeta(1,1,a,b);  
          endif; 
       
            j = j +1;
       endo;
      
      elseif pshape[i] == 2; /* GAMMA PRIOR */
         b = pstdd[i]^2/pmean[i];
         a = pmean[i]/b;
         parasim_t[.,i] = b*rndgam(block_size,1,a);
      
      elseif pshape[i] == 3; /* GAUSSIAN PRIOR */
         a = pmean[i];
         b = pstdd[i];
         parasim_t[.,i] = a + b*rndn(block_size,1);

      elseif pshape[i] == 4; /* INVGAMMA PRIOR */
      
       /*
         a = pmean[i];
         b = pstdd[i];
         parasim_t[.,i] = sqrt( b*a^2/sumc( (rndn(b,block_size))^2 ) );
       */

        j = 1;
        do until j > block_size;
            nn = 2;
            d = pstdd[i]^2/ pi ;
            c = rndc(nn);
            t2 = c/d;
            parasim_t[j,i]= sqrt( 1/t2 );
            j = j +1;
       endo;

      elseif pshape[i] == 5; /* uniform prior */
         a = pmean[i];
         b = pstdd[i];
         parasim_t[.,i] = a + (b-a)*rndu(block_size,1);

      elseif pshape[i] == 0; /* No Prior - Fixed */
         a = pmean[i];
         parasim_t[.,i] = a*ones(block_size,1);

      endif;

      i = i+1;
   endo; 


   locate 1,1;
   "Processing Block" block_it;   @ parasim_t[.,9:18]';@
   
   wr = writer(fhtheta,parasim_t);

   block_it = block_it + 1;

endo;

closeall fhtheta;
end;

#include c:\dge_sw\ceemod.src;


/**************************************/

proc rndc(a);       /* chisquare */
local x,w;

a = a/2;
w= rg1(a);
x = w * 2;

retp(x);
endp;


proc rndb(a,b);                /* beta */
local x,a1n,a1d;

a1n = rg1(a);
a1d = rg1(b);
x = a1n / (a1n + a1d);

retp(x);
endp;


proc rg1(a);
local x,j,w,u;

if a gt 1;
x = rg2(a);
elseif a lt 1;
a = a + 1;
u = rndu(1,1);
x = rg2(a)*u^(1/a);
elseif a == 1;
x = -ln(rndu(1,1));
endif;

retp(x);
endp;


proc rg2(a);
local gam,accept,b,c,j,x,z,u,v,w,y;

b = a-1;
c = 3*a - .75;
accept = 0;
do while accept == 0;
u = rndu(1,1);
v = rndu(1,1);
w = u*(1-u);
y = sqrt(c/w)*(u-.5);
x = b+y;
if x ge 0;
z = 64*(w^3)*(v^2);
accept = z le ( 1-(2*y^2)/x );
if accept == 0;
accept = ln(z) le 2*(b*ln(x/b) - y);
endif;
endif;
endo;

retp(x);
endp;
