function output=finalstateExpGtrellis(val)


%final states (individual
for n = 1:val*2
   
    output(n,2)=n+2*val;
    output(n,3)=1;
    output(n,1)=n;
end

%recursive states
for n = val*2+1:val*4
    output(n,2)=4*val+1;
    output(n,3)=0;
    output(n,1)=n-val*2;
end

output(val*4+1,1)=4*val+1;
output(val*4+1,2)=4*val+2;
output(val*4+1,3)=1;

output(val*4+2,1)=4*val+1;
output(val*4+2,2)=4*val+2;
output(val*4+2,3)=0;

output(val*4+3,1)=4*val+2;
output(val*4+3,2)=4*val+1;
output(val*4+3,3)=0;

output(val*4+4,1)=4*val+2;
output(val*4+4,2)=4*val+3;
output(val*4+4,3)=1;

% output(val*4+5,1)=4*val+3;
% output(val*4+5,2)=4*val+4;
% output(val*4+5,3)=0;
% 
% output(val*4+6,1)=4*val+3;
% output(val*4+6,2)=4*val+5;
% output(val*4+6,3)=1;

