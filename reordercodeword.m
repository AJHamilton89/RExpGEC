function [newcodeword] = reordercodewordRExpG(codeword,k)

% newcodeword = zeros(1,length(codeword));
% codewordflip=fliplr(codeword);
% 
% if length(codeword)<=2
%     newcodeword(1)=codeword(1);
%     newcodeword(2)=codeword(2);
%     
% else
%     newcodeword(1)=codewordflip(1);
%     newcodeword(end-(-1+k))=codewordflip(end); % newcodeword(end-(1-k))=codewordflip(end);
%     for n=2:2:length(codeword)
%         newcodeword(n)=codewordflip(((floor((length(codeword))/2)+(n/2))));%not sure if floor or ceil...
%     end
%     
%     for n=3:2:length(codeword)
%         newcodeword(n)=codewordflip(ceil(n/2));
%     end
%     newcodeword=fliplr(newcodeword);
%     
    %%%
    newcodeword = zeros(1,length(codeword));
    newcodeword(1)=codeword(1);
    if k >=1
        for n=1:k
           newcodeword(end-n+1)=codeword(end-n+1)
        end
    else
    end
    
    transition_section=codeword(2:end-k);
   % length(codeword)-k-1 -- the bit in the middle that needs swapping
    for n=2:2:(length(codeword)-k) %fill up UEC bits
        newcodeword(n)=transition_section(ceil(length(transition_section)/2)+n/2);%codeword(((floor((length(codeword))/2)+(n/2))));%not sure if floor or ceil...
    end
    
    for n=3:2:(length(codeword)-k) % fill up FLC bits
        newcodeword(n)=transition_section(floor(n/2));
    end
   
end
