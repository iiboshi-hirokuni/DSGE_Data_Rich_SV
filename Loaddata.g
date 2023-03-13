/* filename: loaddata.g
** description: reads time series data from ASCII data set into
** Gauss variable SERIES
** created:  12/10/01
** Modifications: 10/12/4 change for BGGGK data(datasel=9(else) ), change data_us.txt for BGGGK data
**
*/

if datasel == 1;

   /* USE SIMULATED DATA
   */

   lsimdata  = datapath $+ "\\" $+ lmodel $+ "sim" $+ dataselstr;
   open fhsimdata = ^lsimdata for read;
   series_YT = readr(fhsimdata,Tsim[datasel]);
   closeall fhsimdata;   
   
   ti = seqa(1,1,Tsim[datasel]);

elseif datasel == 2;  

   /* USE ACTUAL DATA
     */

titlestr =  "Consumption_1" $|"Consumption_2" $|"Consumption_3" $|"Consumption_4" $|
                   "Wage_1" $| "Wage_2" $| "Wage_3" $| "Wage_4" $| 
		           "Labor_1" $| "Labor_2" $| "Labor_3" $| "Labor_4" $|
                   "Nominal_Interest_Rate "  $| 
            	   "Inflation_1" $| "Inflation_2" $| "Inflation_3" $| "Inflation_4" $|  
                   "Output_1" $| "Output_2" $| "Output_3" $| "Output_4" $|
                   "Spread_1"  $|  "Spread_2"  $|  "Spread_3"  $|  "Spread_4"  $|
                   "Bank_Leverage_1" $| "Bank_Leverage_2" $|"Bank_Leverage_3" $| "Bank_Leverage_4" $|
		           "Investment_1" $| "Investment_2" $| "Investment_3" $| "Investment_4" $|
                   "Firm_Leverage_1" $|  "Firm_Leverage_2" $| "Firm_Leverage_3" $|                   
                   "Borrowing_Rate_1" $| "Borrowing_Rate_2" $| "Borrowing_Rate_3" $|
				   "Borrowing_Rate_4"  ; 
 titlestr_save ={};				   

 load path   = ^datapath series_yt[104,41]= "US_data_40.csv";

   ti          = seqa(1985.25,0.25,103);

      series_yt = real(series_yt[2:104,2:41]);
    
   
   if  mhrun == "3" ;
   
       /*
	   series_yt = series_yt[.,1]~series_yt[.,5]~ series_yt[.,9]~ series_yt[.,13]~
	               series_yt[.,14]~series_yt[.,18]~ series_yt[.,22]~ series_yt[.,26]~
	               series_yt[.,30]~series_yt[.,34]~ series_yt[.,37];
       */
       load path = ^priopath case_c_zz_1[40+1,12+1] = "case_a_zz.csv";
 
         case_c_zz_1 = real(case_c_zz_1[2:41,1+1:12+1]);

         series_yt_1 = series_yt;  

         index_type = 1;  case_c_zz={};  series_yt={}; 

         do until index_type > 40;

            if case_c_zz_1[index_type,4] ne 0;
               
               case_c_zz = case_c_zz|case_c_zz_1[index_type,.]; 

               series_yt = series_yt~series_yt_1[.,index_type];
			  
			   titlestr_save = titlestr_save $|titlestr[index_type];	 
 
            endif;

           index_type = index_type + 1;

          endo;              

      
/****************************
 Drawing Graph of Data  
***************************/
                  
   graphset;
    begwind;
    //margin(0,0,0.5,0.5);
    window(4,3,0);           @@
    //fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    //_pcolor = {4 4 4};
    //_pltype = {6 1 1};
    //_protate = 0; 
    //_plctrl = 0;
    //_plwidth = 5.5;
    //_psymsiz = 5;
    //_plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";

       yy = series_yt ;
	   
       y_ind = 1;
       do until y_ind > 11;
       /* iterate over endogenous variables */
          
          title(titlestr_save[y_ind]);

          if (datasel ==2) or (datasel ==3);
              xtics(1985, 2010,2,4);
              x=seqa(0/4,1/4, rows(yy))+1985;
         
          endif;
          
          xy(x, yy[.,y_ind] );  

          nextwind;
                 
          y_ind = y_ind+1;

       endo;
       
    
       endwind;   

  elseif mhrun == "5";

         load path = ^priopath case_c_zz_1[40+1,12+1] = "case_c_zz.csv";
 
         case_c_zz_1 = real(case_c_zz_1[2:41,1+1:12+1]);

         series_yt_1 = series_yt;  

         index_type = 1;  case_c_zz={};  series_yt={}; 

         do until index_type > 40;

            if case_c_zz_1[index_type,4] ne 0;
               
               case_c_zz = case_c_zz|case_c_zz_1[index_type,.]; 

               series_yt = series_yt~series_yt_1[.,index_type];
			  
			   titlestr_save = titlestr_save $|titlestr[index_type]; 
 
            endif;

           index_type = index_type + 1;

          endo;              

    endif;
        
      
endif;

"data setting table" ;
print case_c_zz;  
"press any key!!";
wait;


