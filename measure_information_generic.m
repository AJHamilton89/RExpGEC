function information = measure_information_generic(dict,maxngram)



totalcount=double(sum(dict(:,maxngram+1)));

for n=1:length(dict(:,1))
    p(n)=double(dict(n,maxngram+1))/totalcount;
    entropy(n)=-p(n)*log2(p(n));
end


for n=1:length(dict(:,1))
  numelements(n)=nnz(dict(n,1:maxngram));
end


charpersym=sum(numelements.*p);

information=sum(entropy)/charpersym;