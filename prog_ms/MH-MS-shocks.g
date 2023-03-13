
/*  GIBBS SAMPLING APPROACH TO HAMILTON'S (1989) AR(0) MODEL WITH
    A MARKOV-SWITCHING Vola , APPLIED TO REAL GNP.

    BASED ON THE ORIGINAL CODE WRITTEN BY ALTERT AND CHIB (JBES, 1993).
    -------------------------------------------------------------------

    A MULTI-MOVE GIBBS SAMPLING IS INCORPORATED.
 */


     new;
     #lineson;
     library pgraph;

     format /m1 /rd 9,5;

     data = getpath("bgggk_shocks.csv");
	 print data;
	 
	 load data[103,9]= ^data;
     
	 n = 8;   // number of shocks 
	 /*
	       "Productivity Shock " $|
           "Preference Shock " $|
           "Firm Net Worth Shock" $|
           "Bank Net wortk Shock" $| 
           "Government Spending Shock " $|
           "Investment Shock " $|
           "Labor Supply Shock " $|          
           "Monetary Policy Shock " ;  
	*/	   
									
     ytt = data[2:103,n+1];		//ytt; wait;

     //ytt = ytt - meanc(ytt);

t=rows(ytt);     

   @===== DESIGN  FOR THE GIBBS SAMPLER ===========@

       N0 = 2000;       @NUMBER OF DRAWS TO LEAVE OUT@
       MM = 10000;       @NUMBER OF DRAWS@
       LL = 5;          @EVERY L-TH VALUE IS USED@

       CAPN = N0 + MM;  @TOTAL NUMBER OF DRAWS@


   @===== STORAGE SPACES===========================@

        MU0MM={}; MU1MM={}; SIG20MM={}; SIG21MM={}; PMM={};QMM={};

        sttmm=zeros(t,1);
		
   @=================  PRIOR PARAMETERS ====================@

     U1_01_=2;  U1_00_=8;    @FOR A1TT@
     U1_10_=1;  U1_11_=9;     @FOR B1TT@

     R0_M=EYE(2)*1;     @For MU0, MU1@
     T0_M=0|1;

     D0_=0;
     V0_=0;              @FOR SIG_0^2@		
		
  @========INITIAL VALUES FOR THE GIBBS SAMPLER ===========@

           A1TT=1-0.8;
           B1TT=1-0.8;
		   pq = a1tt|b1tt;
		   
           //MU0TT=0; MU1TT=0;
		   
           SIG20TT=0.5; SIG21TT=0.5;   
		   para = sig20TT|sig21tt;

  ITR=1;
  DO WHILE ITR LE CAPN;      ITR;;"th  iteration.......";

  @=================START SAMPLING================@
  
      {stt, para, pq} = ms_sampler(ytt, para, pq);

       PTT=1-pq[2]; //B1TT;
       QTT=1-pq[1]; //A1TT;

  @===============================================@

                       locate 10,1;
                       //print "this is tranmat" tranmat;
                       //"mu";; mu0tt;;mu1tt;
                       "SIG20 and Sig21";; sig20tt;; SIG21TT;
                       "P AND Q";; PTT;;QTT;

       IF ITR GT N0;
	   
	       SiG20TT= para[1]; SiG21TT= para[2]; 
	   
           SIG20MM=SIG20MM~SIG20TT;
           SIG21MM=SIG21MM~SIG21TT; 
		   
		   //MU0MM=MU0MM~MU0TT; MU1MM=MU1MM~MU1TT; 
		   PMM=PMM~PTT; QMM=QMM~QTT;

           STTMM=STTMM+STT;


       ENDIF;

  ITR=ITR+1;
  ENDO;
  
       STTMM=STTMM/MM;
       FNL_MAT=SIG20MM'~SIG21MM'~MU0MM'~MU1MM'~PMM'~QMM';
       SAVE FNL_MAT,STTMM;

       output file=gibs_ms0.dta reset; sttmm; output off;


     @======== calculation of posterior expectations and std deviations@

          INDX = ZEROS(MM,1);
          INDX[1] = 1;
          J = 1;
          DO WHILE J LE MM;
               INDX[J] = 1;
          J = J + LL;
          ENDO;


          SORT_OUT=ZEROS(cols(fnl_mat),3);
          I_NDX=1;
          DO UNTIL I_NDX>COLS(FNL_MAT);

                tmpm1 = selif(FNL_MAT[.,I_NDX],indx);
                MN_OUT = meanc(tmpm1);
                STD_OUT = stdc(tmpm1);
                MED_OUT = median(tmpm1);
                u2out = MN_OUT~STD_OUT~MED_OUT;

        SORT_OUT[I_NDX,.]=U2OUT;

           I_NDX=I_NDX+1;
       ENDO;

OUTPUT FILE=GIBS_MS0.OUT RESET;

   "FOR SIG2, MU0, MU1, P, Q";
   "mean ,  STD, Median,  "
   "================================================";

    SORT_OUT;
   "PR[S_t=1|I(T)]";
   "================================================";
    STTMM;
	
    //OUTPUT OFF;

    begwind;
    margin(0,0,0.2,0.2);
    window(2,1,0);
	 fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    _pcolor = {2 3 3 3};
    _pltype = {6 1 1 1};
    _protate = 0; 
    _plctrl = 0;
    _plwidth = 5.5;
    _psymsiz = 5;
    _plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";
	
	 xtics(1985, 2010,2,4);
     x=seqa(0/4,1/4, rows(ytt))+1985;
	 
	 titlestr = "Productivity Shock " $|
           "Preference Shock " $|
           "Firm Net Worth Shock" $|
           "Bank Net wortk Shock" $| 
           "Government Spending Shock " $|
           "Investment Shock " $|
           "Labor Supply Shock " $|          
           "Monetary Policy Shock " ; 

title(titlestr[n]);
xy(x,ytt); 

nextwind;

title("Probability of Regime of High Volatility");
xy(x,1-sttmm);

endwind;
	
   
END;

#include c:\dsge\dsge_kk\prog_ms\ms-sampler.src
