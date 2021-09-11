function information = measure_information_ten(dict)



totalcount=double(sum(dict(:,11)));

for n=1:length(dict)
    p(n)=double(dict(n,11))/totalcount;
    entropy(n)=-p(n)*log2(p(n));
end


for n=1:length(dict)
  numelements(n)=nnz(dict(n,1:10));
end


charpersym=sum(numelements.*p);

information=sum(entropy)/charpersym;