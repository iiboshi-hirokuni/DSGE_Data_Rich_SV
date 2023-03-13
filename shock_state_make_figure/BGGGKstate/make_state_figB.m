load state_B.txt
BBB=state_B;

BB=zeros(100,4*21);
for i=1:21
BB(:,1+(i-1)*4:4+(i-1)*4)=BBB((i-1)*100+1:(i-1)*100+100,1:4);
end

figure
subplot(7,3,7)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,4)-i*(4*BB(:,4)-4*BB(:,3)),'c','LineWidth',2)    
%plot(1985.50:0.25:2010.25,BB(:,3),'r--')
%plot(1985.50:0.25:2010.25,BB(:,4),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,1),'black')
plot(1985.50:0.25:2010.25,4*BB(:,2),'b','LineWidth',2)
hold off
title('Nominal Interest Rate')


subplot(7,3,1)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,8)-i*(BB(:,8)-BB(:,7)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,7),'r--')
%plot(1985.50:0.25:2010.25,BB(:,8),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,5),'black')
plot(1985.50:0.25:2010.25,BB(:,6),'b','LineWidth',2)
hold off
title('Output')

subplot(7,3,2)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,12)-i*(BB(:,12)-BB(:,11)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,11),'r--')
%plot(1985.50:0.25:2010.25,BB(:,12),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,9),'black')
plot(1985.50:0.25:2010.25,BB(:,10),'b','LineWidth',2)
hold off
title('Consumption')

subplot(7,3,3)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,16)-i*(BB(:,16)-BB(:,15)),'c','LineWidth',2)       
%plot(1985.50:0.25:2010.25,BB(:,15),'r--')
%plot(1985.50:0.25:2010.25,BB(:,16),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,13),'black')
plot(1985.50:0.25:2010.25,BB(:,14),'b','LineWidth',2)
hold off
title('Investment')

subplot(7,3,4)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,20)-i*(4*BB(:,20)-4*BB(:,19)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,19),'r--')
%plot(1985.50:0.25:2010.25,BB(:,20),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,17),'black')
plot(1985.50:0.25:2010.25,4*BB(:,18),'b','LineWidth',2)
hold off
title('Inflation')

subplot(7,3,5)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,24)-i*(BB(:,24)-BB(:,23)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,23),'r--')
%plot(1985.50:0.25:2010.25,BB(:,24),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,21),'black')
plot(1985.50:0.25:2010.25,BB(:,22),'b','LineWidth',2)
hold off
title('Wage')

subplot(7,3,6)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,28)-i*(BB(:,28)-BB(:,27)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,27),'r--')
%plot(1985.50:0.25:2010.25,BB(:,28),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,25),'black')
plot(1985.50:0.25:2010.25,BB(:,26),'b','LineWidth',2)
hold off
title('Labor')

subplot(7,3,8)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,32)-i*(4*BB(:,32)-4*BB(:,31)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,31),'r--')
%plot(1985.50:0.25:2010.25,BB(:,32),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,29),'black')
plot(1985.50:0.25:2010.25,4*BB(:,30),'b','LineWidth',2)
hold off
title('Corporate Borrowing Rate')

subplot(7,3,11)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,36)-i*(BB(:,36)-BB(:,35)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,35),'r--')
%plot(1985.50:0.25:2010.25,BB(:,36),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,33),'black')
plot(1985.50:0.25:2010.25,BB(:,34),'b','LineWidth',2)
hold off
title('Bank Leverage')

subplot(7,3,10)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,40)-i*(BB(:,40)-BB(:,39)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,39),'r--')
%plot(1985.50:0.25:2010.25,BB(:,40),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,37),'black')
plot(1985.50:0.25:2010.25,BB(:,38),'b','LineWidth',2)
hold off
title('Corporate Leverage')

subplot(7,3,9)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,44)-i*(4*BB(:,44)-4*BB(:,43)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,41),'black')
plot(1985.50:0.25:2010.25,4*BB(:,42),'b','LineWidth',2)
hold off
title('External Premium')

subplot(7,3,12)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,48)-i*(BB(:,48)-BB(:,47)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,45),'black')
plot(1985.50:0.25:2010.25,BB(:,46),'b','LineWidth',2)
hold off
title('Consumption 2')

subplot(7,3,13)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,52)-i*(BB(:,52)-BB(:,51)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,49),'black')
plot(1985.50:0.25:2010.25,BB(:,50),'b','LineWidth',2)
hold off
title('Investment 2')

subplot(7,3,14)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,56)-i*(4*BB(:,56)-4*BB(:,55)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,53),'black')
plot(1985.50:0.25:2010.25,4*BB(:,54),'b','LineWidth',2)
hold off
title('Inflation 2')

subplot(7,3,15)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,60)-i*(4*BB(:,60)-4*BB(:,59)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,57),'black')
plot(1985.50:0.25:2010.25,4*BB(:,58),'b','LineWidth',2)
hold off
title('Inflation 3')

subplot(7,3,16)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,64)-i*(BB(:,64)-BB(:,63)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,61),'black')
plot(1985.50:0.25:2010.25,BB(:,62),'b','LineWidth',2)
hold off
title('Labor 2')

subplot(7,3,17)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,68)-i*(BB(:,68)-BB(:,67)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,65),'black')
plot(1985.50:0.25:2010.25,BB(:,66),'b','LineWidth',2)
hold off
title('Labor 3')

subplot(7,3,20)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,72)-i*(BB(:,72)-BB(:,71)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,69),'black')
plot(1985.50:0.25:2010.25,BB(:,70),'b','LineWidth',2)
hold off
title('Bank Leverage 2')

subplot(7,3,21)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,BB(:,76)-i*(BB(:,76)-BB(:,75)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,BB(:,73),'black')
plot(1985.50:0.25:2010.25,BB(:,74),'b','LineWidth',2)
hold off
title('Bank Leverage 3')

subplot(7,3,18)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,80)-i*(4*BB(:,80)-4*BB(:,79)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,77),'black')
plot(1985.50:0.25:2010.25,4*BB(:,78),'b','LineWidth',2)
hold off
title('External Premium 2')

subplot(7,3,19)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.25,4*BB(:,84)-i*(4*BB(:,84)-4*BB(:,83)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.25,BB(:,43),'r--')
%plot(1985.50:0.25:2010.25,BB(:,44),'r--')
end
plot(1985.50:0.25:2010.25,4*BB(:,81),'black')
plot(1985.50:0.25:2010.25,4*BB(:,82),'b','LineWidth',2)
hold off
title('External Premium 3')

