function output=buildgenericUECtree(k)

for n = 1:2*(2^k-1)
   
    output(n,1)=ceil(n/2);
    output(n,2)=n+1;
    output(n,3)=mod(n,2);
%     if n==1
%         
%     
%     output(n,2)=2*n;
%     output(n,3)=1;
%         
%     elseif mod(n,2)==0
%     
%        
%     output(n,2)=2*n;
%     output(n,3)=1;
%     
%     else
%    
%         output(n,2)=(2*n-1)+1
%         output(n,3)=0;
%     
%     
%     end
end
end