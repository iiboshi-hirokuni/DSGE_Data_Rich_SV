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


titlename={'TFP Shock (ε^A)' ;
           'Preference Shock (ε^C)' ; 
           'Corporate Net Worth Shock (ε^E)' ;
           'Bank Net Worth Shock (ε^F)';            
           'Government Spending Shock (ε^G)' ;                          
           'Investment Specific Technology Shock (ε^K)' ;
           'Labor Supply Shock (ε^L)';
           'Monetary Policy Shock (ε^R)'};

titlename_d={'TFP Shock (ε^A)' ;
           'Preference Shock (ε^C)' ; 
           'Corporate Net Worth Shock (ε^E)' ;
           'Bank Net Worth Shock (ε^F)';            
           'Government Spending Shock (ε^G)' ;                          
           'Investment Specific Technology Shock (ε^K)' ;
           'Labor Supply Shock (ε^L)';
           'Monetary Policy Shock (ε^R)'};     
       
       
legend_name={'Case C'; 'Case D';'Case D'; 'Case D'};      
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end
     
const_vola = ones(T,1)*[0.564;0.238;0.787;0.228]';

figure('Name','SV case','File','SV_case_c.fig',...
       'Position',[100,100,450,800]);
    
    date_end = 1985.5+(T-1)*0.25;
    
    %%  select macro varibales
    control=[1;3;4;8; 2; 6; 7;5];
    
    case_A = 1; 
    case_B = 0;
    
for k=1:8  
   subplot(8,1,k)
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
    
    if k < 5
    plot(1985.50:0.25:date_end,const_vola(:,k),'r--','LineWidth',2)
    end  
    
        %legend(legend_name(1),'Location','NorthWest')
        title( titlename(j),'FontSize',18 )  
  if k == 3
       ylim([0 2])
  elseif k == 1
       ylim([0 1.2])  
  elseif k == 2
       ylim([0 0.5])
  elseif k == 4
       ylim([0 0.3])     
  elseif k == 6
       ylim([0 1.2])  
  elseif k == 7
       ylim([0 3])    
  elseif k == 8
       ylim([0 1.2])     
  end    
                   
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',10);
  hold off
      
end

figure('Name','SV and Data Rich case','File','SV_case_d.fig',...
      'Position',[600,100,450,800])

case_A = 0;
case_B = 1 ;
    
for k=1:8  
   subplot(8,1,k)
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
%           set(h(2),'FaceColor',[1 0.4 0.6])      
          set(h(2),'FaceColor',[0.5 1 1])  
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value  
       
  end
 
   if case_A == 1      
      plot(1985.50:0.25:date_end,AA(:,1+3*(j-1)),'b','LineWidth',2)
     
   end   
   
    if case_B == 1      
     plot(1985.50:0.25:date_end,BB(:,1+3*(j-1)),'b','LineWidth',2)
    end
    
    if k < 5
    plot(1985.50:0.25:date_end,const_vola(:,k),'r--','LineWidth',2)
    end  
    
        %legend(legend_name(1),'Location','NorthWest')
        title( titlename(j),'FontSize',18 )  
        
                   
  xlim([min(1985.25),max(2013.0)])
  if k == 3
       ylim([0 2])
  elseif k == 1
       ylim([0 1.2])  
  elseif k == 2
       ylim([0 0.5])
  elseif k == 4
       ylim([0 0.3])     
  elseif k == 6
       ylim([0 1.2])  
  elseif k == 7
       ylim([0 3])    
  elseif k == 8
       ylim([0 1.2])     
  end    
  set(gca,'fontsize',10);
  hold off
      
end


 

