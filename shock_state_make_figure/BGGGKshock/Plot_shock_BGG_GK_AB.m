case_A = 1 %  on-> 1,  off-> 0
case_B = 0 %  on-> 1,  off-> 0

n = 8;  % number of shocks 
T = 102;  % number of periods

load shock_A.txt
AAA=shock_A;
AA=zeros(T,3*n);

load shock_B.txt
BBB=shock_B;
BB=zeros(T,3*n);


titlename={'TFP Shock (É√A)' ;
           'Preference Shock (É√C)' ; 
           'Entrepreneur Net Worth Shock (É√E)' ;
           'Banking Sector Net Worth Shock (É√F)';            
           'Government Spending Shock (É√G)' ;                          
           'Investment Specific Technology Shock (É√K)' ;
           'Labor Supply Shock (É√L)';
           'Monetary Policy Shock (É√R)'};
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end


    figure(1)
    
for j=1:n  
   subplot(3,3,j)
   hold on
 
   if case_A == 1
      for i=0:0.01:1
         plot(1985.50:0.25:2010.75,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
      end
   
         plot(1985.50:0.25:2010.75,AA(:,1+3*(j-1)),'b','LineWidth',2)
   end
         
         
  if case_B == 1          
     for i=0:0.01:1
        plot(1985.50:0.25:2010.75,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
     end
         
        plot(1985.50:0.25:2010.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
  
        
        title( titlename(j),'FontSize',10 )  
            %title('éYèo','FontSize',14)
  xlim([min(1985.25),max(2011.0)])
  set(gca,'fontsize',12);
  hold off
      
end
 

