load shock_B.txt
BBB=shock_B;
BB=BBB;
BB=zeros(100,3*8);
for i=1:8
BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*100+1:(i-1)*100+100,1:3);
end

%plot(1981.00:0.25:1995.50,BB(:,3)-i*(BB(:,3)-BB(:,2)),'c')

figure
subplot(3,3,1)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,3)-i*(BB(:,3)-BB(:,2)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,1),'b','LineWidth',2)
hold off
title('TFP Shock (εA)')

subplot(3,3,2)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,6)-i*(BB(:,6)-BB(:,5)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,4),'b','LineWidth',2)
hold off
title('Preference Shock (εC)')

subplot(3,3,3)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,9)-i*(BB(:,9)-BB(:,8)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,7),'b','LineWidth',2)
hold off
title('Entrepreneur Net Worth Shock (εE)')

subplot(3,3,4)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,12)-i*(BB(:,12)-BB(:,11)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,10),'b','LineWidth',2)
hold off
title('Banking Sector Net Worth Shock (εF)')

subplot(3,3,5)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,15)-i*(BB(:,15)-BB(:,14)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,13),'b','LineWidth',2)
hold off
title('Government Spending Shock (εG)')

subplot(3,3,6)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,18)-i*(BB(:,18)-BB(:,17)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,16),'b','LineWidth',2)
hold off
title('Investment Specific Technology Shock (εK)')

subplot(3,3,7)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,21)-i*(BB(:,21)-BB(:,20)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,19),'b','LineWidth',2)
hold off
title('Labor Supply Shock (εL)')

subplot(3,3,8)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,24)-i*(BB(:,24)-BB(:,23)),'c','LineWidth',2)
end
plot(1985.50:0.25:2010.25,BB(:,22),'b','LineWidth',2)
hold off
title('Monetary Policy Shock (εR)')


return






%%%%%%1970.00 1981.00   1998.50 1995.75           %CaseA 1995.75?%1998Q3なのでなく、データ部分が1970Q1を取れていない様子?ちがう？いずれにしても一期間ずれている？
subplot(3,3,1)
hold on                                    
plot(1981.00:0.25:1995.75,AA(:,1),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,1),'b') %SW
hold off
title('Preference Shock')

%subplot(3,3,2)
%hold on
%plot(1981.00:0.25:1995.50,AA(:,2),'g--','LineWidth',1)
%plot(1981.00:0.25:1995.75,BB(:,2),'b') %SW
%hold off
%title('Wage Markup Shock')

%subplot(3,3,3)
%hold on
%plot(1981.00:0.25:1995.50,AA(:,3),'g--','LineWidth',1)
%plot(1981.00:0.25:1995.75,BB(:,3),'b') %SW
%hold off
%title('Price Markup Shock')

subplot(3,3,2)
hold on
plot(1981.00:0.25:1995.75,AA(:,2),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,4),'b') %SW
hold off
title('Investment Shock')

subplot(3,3,3)
hold on
plot(1981.00:0.25:1995.75,AA(:,3),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,5),'b') %SW
hold off
title('Labor Supply Shock')

subplot(3,3,4)
hold on
plot(1981.00:0.25:1995.75,AA(:,4),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,6),'b') %SW
hold off
title('Wage Markup Shock')

subplot(3,3,5)
hold on
plot(1981.00:0.25:1995.75,AA(:,5),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,7),'b') %SW
hold off
title('Productivity Shock')

subplot(3,3,6)
hold on
plot(1981.00:0.25:1995.75,AA(:,6),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,8),'b') %SW
hold off
title('Price Markup Shock')

subplot(3,3,7)
hold on
plot(1981.00:0.25:1995.75,AA(:,7),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,9),'b') %SW
hold off
title('Government Spending Shock')

subplot(3,3,8)
hold on
plot(1981.00:0.25:1995.75,AA(:,8),'r','LineWidth',2)
%plot(1981.00:0.25:1995.75,BB(:,9),'b') %SW
hold off
title('Monetary Policy Shock ')