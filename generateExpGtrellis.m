function output=generateExpGtrellis(val)

for n=1:val*2
    if n<=val
        output(n,2)=1; % switcheroo
        output(n,1)=-val;
    else
        output(n,2)=0;
        output(n,1)=-2*val;
    end
    
end

for n=val*2+1:val*4
    if mod(n,2) ==1
        output(n,2)=0; 
    else
        output(n,2)=1;
    end
    output(n,1)=-ceil((n-1)/2);
%     output
end

% for
end
        