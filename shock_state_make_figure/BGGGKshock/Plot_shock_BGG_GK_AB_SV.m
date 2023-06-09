%  �O���t�̏o�͐ݒ�

case_A = 1 %  print on-> 1,  print off-> 0
case_B = 1 %  print on-> 1,  print off-> 0

% �ʃV���b�N�̃O���t�@Figure 2�@�̐ݒ�
shock_no = 8;  % 4:  'Banking Sector Net Worth Shock (��F)';  
               % 3:   'Entrepreneur Net Worth Shock (��E)' ;
               % 8:    'Monetary Policy Shock (��R)'};
               
%****************************************************************               

n = 8;  % number of shocks 
T = 102;  % number of periods

load sv_vola_sample_A.txt
AAA=sv_vola_sample_A;  
AA=zeros(T,3*n);

load sv_vola_sample_C.txt
BBB=sv_vola_sample_C;
BB=zeros(T,3*n);


titlename={'TFP Shock (��A)' ;
           'Preference Shock (��C)' ; 
           'Entrepreneur Net Worth Shock (��E)' ;
           'Banking Sector Net Worth Shock (��F)';            
           'Government Spending Shock (��G)' ;                          
           'Investment Specific Technology Shock (��K)' ;
           'Labor Supply Shock (��L)';
           'Monetary Policy Shock (��R)'};
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end


    figure(1)
    
for j=1:n  
   subplot(3,3,j)
   hold on
 
  
         
         
  if case_B == 1          
     for i=0:0.01:1
        plot(1985.50:0.25:2010.75,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
     end
         
        plot(1985.50:0.25:2010.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
  
  
   if case_A == 1
      for i=0:0.01:1
         plot(1985.50:0.25:2010.75,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
      end
   
         plot(1985.50:0.25:2010.75,AA(:,1+3*(j-1)),'b','LineWidth',2)
   end
  
   if case_B == 1      
     plot(1985.50:0.25:2010.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
   end
   
        title( titlename(j),'FontSize',10 )  
            %title('�Y�o','FontSize',14)
  xlim([min(1985.25),max(2011.0)])
  set(gca,'fontsize',12);
  hold off
      
end

 figure(2)
 
 j = shock_no;
 
  subplot(1,1,1)
   hold on
 
           
  if case_B == 1          
     for i=0:0.01:1
        plot(1985.50:0.25:2010.75,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
     end
         
        plot(1985.50:0.25:2010.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
  
  
   if case_A == 1
      for i=0:0.01:1
         plot(1985.50:0.25:2010.75,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
      end
   
         plot(1985.50:0.25:2010.75,AA(:,1+3*(j-1)),'b','LineWidth',2)
   end
  
   if case_B == 1      
     plot(1985.50:0.25:2010.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
   end
   
         
    
        title( titlename(j),'FontSize',10 )  
            %title('�Y�o','FontSize',14)
  xlim([min(1985.25),max(2011.0)])
  set(gca,'fontsize',12);
  hold off

 

