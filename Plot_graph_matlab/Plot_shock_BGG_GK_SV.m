%  グラフの出力設定

case_A = 1 %  print on-> 1,  print off-> 0
case_B = 1 %  print on-> 1,  print off-> 0

% 個別ショックのグラフ　Figure 2　の設定
shock_no = 8;  % 4:  'Banking Sector Net Worth Shock (εF)';  
               % 3:   'Entrepreneur Net Worth Shock (εE)' ;
               % 8:    'Monetary Policy Shock (εR)'};
               
%****************************************************************               

n = 8;  % number of shocks 
T = 108;  % number of periods

load sv_vola_sample_C.txt
AAA=sv_vola_sample_C;  
AA=zeros(T,3*n);

load sv_vola_sample_D.txt
BBB=sv_vola_sample_D;
BB=zeros(T,3*n);


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
end
     
const_vola = ones(T,1)*[0.564;0.238;0.787;0.228]';

figure(5)
    
    date_end = 1985.5+(T-1)*0.25;
    
    control=[1;3;4;8];
    
for k=1:4  
   subplot(2,2,k)
   hold on  
   
   j=control(k);
         
   %j=k;
      
   if case_A == 1
%       for i=0:0.01:1
%          plot(1985.50:0.25:date_end,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end
      h=area(1985.50:0.25:date_end,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
   
       
   end
   
    if case_B == 1          
%      for i=0:0.01:1
%         plot(1985.50:0.25:date_end,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%      end     

       h=area(1985.50:0.25:date_end,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value  
       
  end
 
   if case_A == 1      
      plot(1985.50:0.25:date_end,AA(:,1+3*(j-1)),'b','LineWidth',2)
     
   end   
   
    if case_B == 1      
     plot(1985.50:0.25:date_end,BB(:,1+3*(j-1)),'r','LineWidth',2)
    end
    
    
    plot(1985.50:0.25:date_end,const_vola(:,k),'--k','LineWidth',2)
      
    
        %legend(legend_name(1),'Location','NorthWest')
        title( titlename(j),'FontSize',14 )  
        
                   
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',12);
  hold off
      
end

 

