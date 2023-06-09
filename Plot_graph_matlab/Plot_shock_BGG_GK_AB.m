%  �O���t�̏o�͐ݒ�

case_A = 1 %  print on-> 1,  print off-> 0
case_B = 1 %  print on-> 1,  print off-> 0

% �ʃV���b�N�̃O���t�@Figure 2�@�̐ݒ�
shock_no = 8;  % 4:  'Banking Sector Net Worth Shock (��F)';  
               % 3:   'Entrepreneur Net Worth Shock (��E)' ;
               % 8:    'Monetary Policy Shock (��R)'};
               
%****************************************************************

n = 8;  % number of shocks 
T = 108;  % number of periods

load shock_sample_A.txt
AAA=shock_sample_A;
AA=zeros(T,3*n);

load shock_sample_B.txt
BBB=shock_sample_B;
BB=zeros(T,3*n);


titlename={'TFP Shock (��A)' ;
           'Preference Shock (��C)' ; 
           'Corporate Net Worth Shock (��E)' ;
           'Bank Net Worth Shock (��F)';            
           'Government Spending Shock (��G)' ;                          
           'Investment Specific Technology Shock (��K)' ;
           'Labor Supply Shock (��L)';
           'Monetary Policy Shock (��R)'};
       
legend_name={'Case A'; 'Case B';'Case B'; 'Case B'};       
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end


    figure(1)
    
    date_end = 1985.5+(T-1)*0.25;
    
for j=1:n  
   subplot(3,3,j)
   hold on
 
   if case_A == 1
%       for i=0:0.01:1          
%           plot(1985.50:0.25:date_end,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end

      h=area(1985.50:0.25:date_end,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
   
   end
   
     if case_A == 1      
      plot(1985.50:0.25:date_end,AA(:,1+3*(j-1)),'b','LineWidth',2)
     
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
 
   
    if case_B == 1      
     plot(1985.50:0.25:date_end,BB(:,1+3*(j-1)),'r','LineWidth',2)
   end
  
  
       % legend(legend_name(1),legend_name(3))
%        legend(legend_name(1),'Location','NorthWest')
        title( titlename(j),'FontSize',12 )  
            %title('�Y�o','FontSize',14)
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',12);
  hold off
      
end
 
  figure(2)
    
for j=1:n  
   subplot(3,3,j)
   hold on
   
    if case_B == 1          
%      for i=0:0.01:1
%         plot(1985.50:0.25:date_end,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%      end

         h=area(1985.50:0.25:date_end,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
         
    end
  
    
      if case_B == 1      
     plot(1985.50:0.25:date_end,BB(:,1+3*(j-1)),'r','LineWidth',2)
   end
   
 
   if case_A == 1
%       for i=0:0.01:1          
%           plot(1985.50:0.25:date_end,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end

       h=area(1985.50:0.25:date_end,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value         
   end  
   
 
   if case_A == 1      
      plot(1985.50:0.25:date_end,AA(:,1+3*(j-1)),'b','LineWidth',2)
     
   end  
          
       % legend(legend_name(1),legend_name(3))
%        legend(legend_name(2),'Location','NorthWest')
        title( titlename(j),'FontSize',12 )  
            %title('�Y�o','FontSize',14)
  xlim([min(1985.25),max(2013.0)])
  set(gca,'fontsize',12);
  hold off
      
end
