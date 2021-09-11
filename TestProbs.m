val1=0;
val2=0;
clear overall1
clear overall2
figure (2)
clf('reset') 

for n = 1:1000000000
    
x=rand;

if rand>0.8
   
    val1=val1+1;
    
else
    
    val2=val2+1;
    
end

overall1(n)=((val1/n)-0.2)^2;
overall2(n)=((val2/n)-0.8)^2;

if overall1(n)==0.8
    pause
end
% 
% printtest=n/100;
% 
% if rem(printtest,1) == 0
%             
%         
% hold on
% semilogx(overall1);
% drawnow
% semilogx(overall2);
% drawnow
% end

end

hold off
% 
loglog(overall1);
hold
loglog(overall2);
