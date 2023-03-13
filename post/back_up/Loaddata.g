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
      load path   = ^datapath series_yt[103,11]= "newdata_us.txt";//20110407series_yt[101,11]
      ti          = seqa(1985.25,0.25,103);
      series_yt = series_yt[1:103,1:11];
   /* USE ACTUAL DATA
     */

   /*load path   = ^datapath series_yt[103,12]= "US_data.txt";

   ti          = seqa(1985.25,0.25,101);

      series_yt = series_yt[3:103,1:11];*/
    
   /*     
   if subT == 2;
      series_yt = series_yt[1:78,.];
      ti        = ti[1:78,.];
   elseif subT == 3;
      series_yt = series_yt[79:152,.];
      ti        = ti[79:152,.];
   elseif subT == 4;
      series_yt = series_yt[92:152,.];
      ti        = ti[92:152,.];
   endif;
   */
      
endif;




