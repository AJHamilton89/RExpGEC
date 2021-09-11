%Test 2nd FLC transition

max = 100;

for l=1:max % need to ignore first two as these are never met
lFLC(l)=ceil((log2(l+2^k)));
FLC(l)=((l+2^k));
end

 A=FLC-(2.^(lFLC-1))-(2.^(lFLC-2));
 
 A>0