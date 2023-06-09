%  OtĢoĶŻč

case_A = 1 %  print on-> 1,  print off-> 0
case_B = 1 %  print on-> 1,  print off-> 0

% ĀŹVbNĢOt@Figure 2@ĢŻč
shock_no = 8;  % 4:  'Banking Sector Net Worth Shock (ĆF)';  
               % 3:   'Entrepreneur Net Worth Shock (ĆE)' ;
               % 8:    'Monetary Policy Shock (ĆR)'};
               
%****************************************************************               

n = 8;  % number of shocks 
T = 108;  % number of periods

load shock_sample_B.txt
AAA=shock_sample_B;   
AA=zeros(T,3*n);

load shock_sample_D.txt
BBB=shock_sample_D;  
BB=zeros(T,3*n);


titlename={'TFP Shock (Ć^A)' ;
           'Preference Shock (Ć^C)' ; 
           'Corporate Net Worth Shock (Ć^E)' ;
           'Bank Net Worth Shock (Ć^F)';            
           'Government Spending Shock (Ć^G)' ;                          
           'Investment Specific Technology Shock (Ć^K)' ;
           'Labor Supply Shock (Ć^L)';
           'Monetary Policy Shock (Ć^R)'};

titlename_d={'TFP Shock (Ć^A)' ;
           'Preference Shock (Ć^C)' ; 
           'Corporate Net Worth Shock (Ć^E)' ;
           'Bank Net Worth Shock (Ć^F)';            
           'Government Spending Shock (Ć^G)' ;                          
           'Investment Specific Technology Shock (Ć^K)' ;
           'Labor Supply Shock (Ć^L)';
           'Monetary Policy Shock (Ć^R)'};     
       
       
legend_name={'Case C'; 'Case D';'Case D'; 'Case D'};      
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end
     
const_vola = ones(T,1)*[0.564;0.238;0.787;0.228]';

figure('Name','CV case','File','Shock_CV.fig',...
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
    
%     if k < 5
%     plot(1985.50:0.25:date_end,const_vola(:,k),'r--','LineWidth',2)
%     end  
    
        %legend(legend_name(1),'Location','NorthWest')
        title( titlename(j),'FontSize',16 )  
  if k == 3
       ylim([-4 2])
  elseif k == 1
       ylim([-3 1])  
  elseif k == 2
       ylim([-2 0.5])
  elseif k == 4
       ylim([-1.5 0.5])  
  elseif k == 5
       ylim([-4 2])       
  elseif k == 6
       ylim([-2 2])  
  elseif k == 7
       ylim([-5 8])    
  elseif k == 8
       ylim([-2 2])     
  end    
                   
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',10);
  hold off
      
end

figure('Name','SV case','File','shock_SV_case.fig',...
       'Position',[600,100,450,800]);

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
    
%     if k < 5
%     plot(1985.50:0.25:date_end,const_vola(:,k),'r--','LineWidth',2)
%     end  
    
        %legend(legend_name(1),'Location','NorthWest')
        title( titlename(j),'FontSize',16 )  
        
                   
  xlim([min(1985.25),max(2013.0)])
 if k == 3
       ylim([-4 2])
  elseif k == 1
       ylim([-3 1])  
  elseif k == 2
       ylim([-2 0.5])
  elseif k == 4
       ylim([-1.5 0.5])  
  elseif k == 5
       ylim([-4 2])       
  elseif k == 6
       ylim([-2 2])  
  elseif k == 7
       ylim([-5 8])    
  elseif k == 8
       ylim([-2 2])     
  end    
 
 set(gca,'fontsize',10);
  hold off
      
end


 

