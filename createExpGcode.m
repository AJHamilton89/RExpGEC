function codeword=createExpGcode(l,k)

codeword = [zeros(1,length(de2bi(l+2^k-1))-1-k),fliplr(de2bi(l+2^k-1))];