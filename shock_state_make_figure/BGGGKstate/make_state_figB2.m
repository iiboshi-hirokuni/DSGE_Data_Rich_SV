load state_B.txt
BBB=state_B;

BB=zeros(102,4*36);
for i=1:36
BB(:,1+(i-1)*4:4+(i-1)*4)=BBB((i-1)*102+1:(i-1)*102+102,1:4);
end

figure

subplot(4,3,1)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,68)-i*(BB(:,68)-BB(:,67)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,65),'black')
plot(1985.50:0.25:2010.75,BB(:,66),'b','LineWidth',2)
hold off
title('Output 1')

subplot(4,3,2)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,72)-i*(BB(:,72)-BB(:,71)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,69),'black')
plot(1985.50:0.25:2010.75,BB(:,70),'b','LineWidth',2)
hold off
title('Output 2')

subplot(4,3,3)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,76)-i*(BB(:,76)-BB(:,75)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,73),'black')
plot(1985.50:0.25:2010.75,BB(:,74),'b','LineWidth',2)
hold off
title('Output 3')

subplot(4,3,4)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,80)-i*(4*BB(:,80)-4*BB(:,79)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,77),'black')
plot(1985.50:0.25:2010.75,4*BB(:,78),'b','LineWidth',2)
hold off
title('Output 4')

subplot(4,3,5)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,4)-i*(4*BB(:,4)-4*BB(:,3)),'c','LineWidth',2)    
%plot(1985.50:0.25:2010.75,BB(:,3),'r--')
%plot(1985.50:0.25:2010.75,BB(:,4),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,1),'black')
plot(1985.50:0.25:2010.75,4*BB(:,2),'b','LineWidth',2)
hold off
title('Consumption 1')

subplot(4,3,6)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,8)-i*(BB(:,8)-BB(:,7)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,7),'r--')
%plot(1985.50:0.25:2010.75,BB(:,8),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,5),'black')
plot(1985.50:0.25:2010.75,BB(:,6),'b','LineWidth',2)
hold off
title('Consumption 2')

subplot(4,3,7)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,12)-i*(BB(:,12)-BB(:,11)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,11),'r--')
%plot(1985.50:0.25:2010.75,BB(:,12),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,9),'black')
plot(1985.50:0.25:2010.75,BB(:,10),'b','LineWidth',2)
hold off
title('Consumption 3')

subplot(4,3,8)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,16)-i*(BB(:,16)-BB(:,15)),'c','LineWidth',2)       
%plot(1985.50:0.25:2010.75,BB(:,15),'r--')
%plot(1985.50:0.25:2010.75,BB(:,16),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,13),'black')
plot(1985.50:0.25:2010.75,BB(:,14),'b','LineWidth',2)
hold off
title('Consumption 4')

subplot(4,3,9)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,112)-i*(4*BB(:,112)-4*BB(:,111)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,109),'black')
plot(1985.50:0.25:2010.75,4*BB(:,110),'b','LineWidth',2)
hold off
title('Investment 1')

subplot(4,3,10)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,116)-i*(4*BB(:,116)-4*BB(:,115)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,113),'black')
plot(1985.50:0.25:2010.75,4*BB(:,114),'b','LineWidth',2)
hold off
title('Investment 2')

subplot(4,3,11)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,120)-i*(4*BB(:,120)-4*BB(:,119)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,117),'black')
plot(1985.50:0.25:2010.75,4*BB(:,118),'b','LineWidth',2)
hold off
title('Investment 3')

subplot(4,3,12)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,52)-i*(BB(:,52)-BB(:,51)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,49),'black')
plot(1985.50:0.25:2010.75,BB(:,50),'b','LineWidth',2)
hold off
title('Inflation 1')

figure

subplot(4,3,1)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,56)-i*(4*BB(:,56)-4*BB(:,55)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,53),'black')
plot(1985.50:0.25:2010.75,4*BB(:,54),'b','LineWidth',2)
hold off
title('Inflation 2')

subplot(4,3,2)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,60)-i*(4*BB(:,60)-4*BB(:,59)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,57),'black')
plot(1985.50:0.25:2010.75,4*BB(:,58),'b','LineWidth',2)
hold off
title('Inflation 3')

subplot(4,3,3)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,64)-i*(BB(:,64)-BB(:,63)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,61),'black')
plot(1985.50:0.25:2010.75,BB(:,62),'b','LineWidth',2)
hold off
title('Inflation 4')

subplot(4,3,4)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,20)-i*(4*BB(:,20)-4*BB(:,19)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,19),'r--')
%plot(1985.50:0.25:2010.75,BB(:,20),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,17),'black')
plot(1985.50:0.25:2010.75,4*BB(:,18),'b','LineWidth',2)
hold off
title('Wage 1')

subplot(4,3,5)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,24)-i*(BB(:,24)-BB(:,23)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,23),'r--')
%plot(1985.50:0.25:2010.75,BB(:,24),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,21),'black')
plot(1985.50:0.25:2010.75,BB(:,22),'b','LineWidth',2)
hold off
title('Wage 2')

subplot(4,3,6)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,28)-i*(BB(:,28)-BB(:,27)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,27),'r--')
%plot(1985.50:0.25:2010.75,BB(:,28),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,25),'black')
plot(1985.50:0.25:2010.75,BB(:,26),'b','LineWidth',2)
hold off
title('Wage 3')

subplot(4,3,7)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,32)-i*(4*BB(:,32)-4*BB(:,31)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,31),'r--')
%plot(1985.50:0.25:2010.75,BB(:,32),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,29),'black')
plot(1985.50:0.25:2010.75,4*BB(:,30),'b','LineWidth',2)
hold off
title('Wage 4')

subplot(4,3,8)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,36)-i*(BB(:,36)-BB(:,35)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,35),'r--')
%plot(1985.50:0.25:2010.75,BB(:,36),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,33),'black')
plot(1985.50:0.25:2010.75,BB(:,34),'b','LineWidth',2)
hold off
title('Labor 1')

subplot(4,3,9)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,40)-i*(BB(:,40)-BB(:,39)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,39),'r--')
%plot(1985.50:0.25:2010.75,BB(:,40),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,37),'black')
plot(1985.50:0.25:2010.75,BB(:,38),'b','LineWidth',2)
hold off
title('Labor 3')

subplot(4,3,10)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,44)-i*(4*BB(:,44)-4*BB(:,43)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,41),'black')
plot(1985.50:0.25:2010.75,4*BB(:,42),'b','LineWidth',2)
hold off
title('Labor 4')

subplot(4,3,11)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,48)-i*(BB(:,48)-BB(:,47)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,45),'black')
plot(1985.50:0.25:2010.75,BB(:,46),'b','LineWidth',2)
hold off
title('Nominal Interest Rate')


subplot(4,3,12)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,132)-i*(BB(:,132)-BB(:,131)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,129),'black')
plot(1985.50:0.25:2010.75,BB(:,130),'b','LineWidth',2)
hold off
title('Corporate Borrowing Rate 1')

figure

subplot(4,3,1)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,136)-i*(BB(:,136)-BB(:,135)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,133),'black')
plot(1985.50:0.25:2010.75,BB(:,134),'b','LineWidth',2)
hold off
title('Corporate Borrowing Rate 2')

subplot(4,3,2)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,BB(:,140)-i*(BB(:,140)-BB(:,139)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,BB(:,137),'black')
plot(1985.50:0.25:2010.75,BB(:,138),'b','LineWidth',2)
hold off
title('Corporate Borrowing Rate 3')

subplot(4,3,3)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,144)-i*(4*BB(:,144)-4*BB(:,143)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,141),'black')
plot(1985.50:0.25:2010.75,4*BB(:,142),'b','LineWidth',2)
hold off
title('Corporate Borrowing Rate 4')

subplot(4,3,4)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,84)-i*(4*BB(:,84)-4*BB(:,83)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,81),'black')
plot(1985.50:0.25:2010.75,4*BB(:,82),'b','LineWidth',2)
hold off
title('External Premium 1')

subplot(4,3,5)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,88)-i*(4*BB(:,88)-4*BB(:,87)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,85),'black')
plot(1985.50:0.25:2010.75,4*BB(:,86),'b','LineWidth',2)
hold off
title('External Premium 2')

subplot(4,3,6)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,92)-i*(4*BB(:,92)-4*BB(:,91)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,89),'black')
plot(1985.50:0.25:2010.75,4*BB(:,90),'b','LineWidth',2)
hold off
title('External Premium 3')

subplot(4,3,7)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,96)-i*(4*BB(:,96)-4*BB(:,95)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,93),'black')
plot(1985.50:0.25:2010.75,4*BB(:,94),'b','LineWidth',2)
hold off
title('External Premium 4')

subplot(4,3,8)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,124)-i*(4*BB(:,124)-4*BB(:,123)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,121),'black')
plot(1985.50:0.25:2010.75,4*BB(:,122),'b','LineWidth',2)
hold off
title('Corporate Leverage 1')

subplot(4,3,9)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,128)-i*(4*BB(:,128)-4*BB(:,127)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,125),'black')
plot(1985.50:0.25:2010.75,4*BB(:,126),'b','LineWidth',2)
hold off
title('Corporate Leverage 3')

subplot(4,3,10)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,100)-i*(4*BB(:,100)-4*BB(:,99)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,97),'black')
plot(1985.50:0.25:2010.75,4*BB(:,98),'b','LineWidth',2)
hold off
title('Bank Leverage 1')

subplot(4,3,11)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,104)-i*(4*BB(:,104)-4*BB(:,103)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,101),'black')
plot(1985.50:0.25:2010.75,4*BB(:,102),'b','LineWidth',2)
hold off
title('Bank Leverage 2')

subplot(4,3,12)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,4*BB(:,108)-i*(4*BB(:,108)-4*BB(:,107)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,BB(:,43),'r--')
%plot(1985.50:0.25:2010.75,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.75,4*BB(:,105),'black')
plot(1985.50:0.25:2010.75,4*BB(:,106),'b','LineWidth',2)
hold off
title('Bank Leverage 3')














