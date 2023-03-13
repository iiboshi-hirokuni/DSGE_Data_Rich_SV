%  グラフの出力設定

case_A = 1 %  print on-> 1,  print off-> 0
case_B = 1 %  print on-> 1,  print off-> 0
case_C = 1 %  print on-> 1,  print off-> 0
case_D = 1 %  print on-> 1,  print off-> 0

% 個別ショックのグラフ　Figure 2　の設定
shock_no = 4;  % 1:  TFP Shock
               % 4:  'Banking Sector Net Worth Shock (εF)';  
               % 3:   'Entrepreneur Net Worth Shock (εE)' ;
               % 8:    'Monetary Policy Shock (εR)'};
               
%****************************************************************               

n = 8;  % number of shocks 
T = 108;  % number of periods

load shock_sample_A.txt
AAA=shock_sample_A;  
AA=zeros(T,3*n);

load shock_sample_B.txt
BBB=shock_sample_B;
BB=zeros(T,3*n);

load shock_sample_C.txt
CCC=shock_sample_C;  
CC=zeros(T,3*n);

load shock_sample_D.txt
DDD=shock_sample_D;
DD=zeros(T,3*n);


titlename={'TFP Shock (εA)' ;
           'Preference Shock (εC)' ; 
           'Corporate Net Worth Shock (εE)' ;
           'Bank Net Worth Shock (εF)';            
           'Government Spending Shock (εG)' ;                          
           'Investment Specific Technology Shock (εK)' ;
           'Labor Supply Shock (εL)';
           'Monetary Policy Shock (εR)'};
       
legend_name={'Case C'; 'Case D';'Case D'; 'Case D'};      
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
    
    CC(:,1+(i-1)*3:3+(i-1)*3)=CCC((i-1)*T+1:(i-1)*T+T,1:3);

    DD(:,1+(i-1)*3:3+(i-1)*3)=DDD((i-1)*T+1:(i-1)*T+T,1:3);
        
end

     
    date_end = 1985.5+(T-1)*0.25;
    

 figure(7)
 
 subplot(2,2,1)
   
 
 j = shock_no; 
 
   if case_A == 1
%       for i=0:0.01:1
%          plot(1985.50:0.25:date_end,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end

%          area(1985.50:0.25:date_end,AA(:,3+3*(j-1)),'c','LineWidth',2)  
   hold on

      h=area(1985.50:0.25:date_end,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
   
           plot(1985.50:0.25:date_end,AA(:,1+3*(j-1)),'b','LineWidth',2)
          
         
   end
                   
        title( 'CaseA','FontSize',18 )  
        
  %ylim([-1.5,0.5])          
  xlim([min(1985.25),max(2013.0)])
%   set(gca,'fontsize',12);
  hold off
  
  subplot(2,2,2)
   hold on
  
if case_B == 1          
%      for i=0:0.01:1
%         plot(1985.50:0.25:date_end,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%      end

  h=area(1985.50:0.25:date_end,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value

         
        plot(1985.50:0.25:date_end,BB(:,1+3*(j-1)),'r','LineWidth',2)
end   
           title( 'Case B','FontSize',18 )  
            
  xlim([min(1985.25),max(2013.0)])
  %ylim([-1.5,0.5])
  set(gca,'fontsize',12);
  hold off
  
  
 subplot(2,2,3)
   hold on
  
if case_C == 1          
%      for i=0:0.01:1
%         plot(1985.50:0.25:date_end,CC(:,3+3*(j-1))-i*(CC(:,3+3*(j-1))-CC(:,2+3*(j-1))),'c', 'LineWidth',2)    
%      end

          h=area(1985.50:0.25:date_end,[ CC(:,2+3*(j-1))  (CC(:,3+3*(j-1))-CC(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
         
        plot(1985.50:0.25:date_end,CC(:,1+3*(j-1)),'b','LineWidth',2)
end   
   
           title( 'Case C','FontSize',18 ) 
           
            
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',12);
  hold off
  
  
  subplot(2,2,4)
   hold on
   
  if case_D == 1          
%      for i=0:0.01:1
%         plot(1985.50:0.25:date_end,DD(:,3+3*(j-1))-i*(DD(:,3+3*(j-1))-DD(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%      end

          h=area(1985.50:0.25:date_end,[ DD(:,2+3*(j-1))  (DD(:,3+3*(j-1))-DD(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
         
        plot(1985.50:0.25:date_end,DD(:,1+3*(j-1)),'r','LineWidth',2)
  end
     
           title( 'Case D','FontSize',18 )  
            
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',12);
  hold off
  
  

