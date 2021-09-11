function information = measure_information(dict)



totalcount=double(sum(dict(:,5)));

for n=1:length(dict)
    p(n)=double(dict(n,5))/totalcount;
    entropy(n)=-p(n)*log2(p(n));
end


for n=1:length(dict)
  numelements(n)=nnz(dict(n,1:4));
end


charpersym=sum(numelements.*p);

information=sum(entropy)/charpersym;