load state_A.txt
AAA=state_A;

AA=zeros(102,4*8);
for i=1:11
AA(:,1+(i-1)*4:4+(i-1)*4)=AAA((i-1)*102+1:(i-1)*102+102,1:4);
end

figure
subplot(4,3,1)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,24)-i*(AA(:,24)-AA(:,23)),'c','LineWidth',2)    
%plot(1985.50:0.25:2010.75,AA(:,3),'r--')
%plot(1985.50:0.25:2010.75,AA(:,4),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,21),'black')
plot(1985.50:0.25:2010.75,AA(:,22),'b','LineWidth',2)
hold off
title('Output')

subplot(4,3,2)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,4)-i*(AA(:,4)-AA(:,3)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,7),'r--')
%plot(1985.50:0.25:2010.75,AA(:,8),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,1),'black')
plot(1985.50:0.25:2010.75,AA(:,2),'b','LineWidth',2)
hold off
title('Consumption')

subplot(4,3,3)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,36)-i*(AA(:,36)-AA(:,35)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,11),'r--')
%plot(1985.50:0.25:2010.75,AA(:,12),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,33),'black')
plot(1985.50:0.25:2010.75,AA(:,34),'b','LineWidth',2)
hold off
title('Investment')

subplot(4,3,4)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,20)-i*(AA(:,20)-AA(:,19)),'c','LineWidth',2)       
%plot(1985.50:0.25:2010.75,AA(:,15),'r--')
%plot(1985.50:0.25:2010.75,AA(:,16),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,17),'black')
plot(1985.50:0.25:2010.75,AA(:,18),'b','LineWidth',2)
hold off
title('Inflation')

subplot(4,3,5)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,8)-i*(AA(:,8)-AA(:,7)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,19),'r--')
%plot(1985.50:0.25:2010.75,AA(:,20),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,5),'black')
plot(1985.50:0.25:2010.75,AA(:,6),'b','LineWidth',2)
hold off
title('Wage')

subplot(4,3,6)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,12)-i*(AA(:,12)-AA(:,11)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,23),'r--')
%plot(1985.50:0.25:2010.75,AA(:,24),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,9),'black')
plot(1985.50:0.25:2010.75,AA(:,10),'b','LineWidth',2)
hold off
title('Labor')

subplot(4,3,7)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,16)-i*(AA(:,16)-AA(:,15)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,27),'r--')
%plot(1985.50:0.25:2010.75,AA(:,28),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,13),'black')
plot(1985.50:0.25:2010.75,AA(:,14),'b','LineWidth',2)
hold off
title('Nominal Interest Rate')

subplot(4,3,8)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,44)-i*(AA(:,44)-AA(:,43)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,31),'r--')
%plot(1985.50:0.25:2010.75,AA(:,32),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,41),'black')
plot(1985.50:0.25:2010.75,AA(:,42),'b','LineWidth',2)
hold off
title('Corporate Borrowing Rate')

subplot(4,3,9)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,28)-i*(AA(:,28)-AA(:,27)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,35),'r--')
%plot(1985.50:0.25:2010.75,AA(:,36),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,25),'black')
plot(1985.50:0.25:2010.75,AA(:,26),'b','LineWidth',2)
hold off
title('External Premium')

subplot(4,3,10)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,40)-i*(AA(:,40)-AA(:,39)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,39),'r--')
%plot(1985.50:0.25:2010.75,AA(:,40),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,37),'black')
plot(1985.50:0.25:2010.75,AA(:,38),'b','LineWidth',2)
hold off
title('Corporate Leverage')

subplot(4,3,11)
hold on
for i=0:0.01:1
plot(1985.50:0.25:2010.75,AA(:,32)-i*(AA(:,32)-AA(:,31)),'c','LineWidth',2)        
%plot(1985.50:0.25:2010.75,AA(:,43),'r--')
%plot(1985.50:0.25:2010.75,AA(:,44),'r--')
end
plot(1985.50:0.25:2010.75,AA(:,29),'black')
plot(1985.50:0.25:2010.75,AA(:,30),'b','LineWidth',2)
hold off
title('Bank Leverage')

