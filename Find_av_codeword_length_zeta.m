function [averagecodewordlength,entropy] = Find_av_codeword_length_zeta(maxcodes,k,s)

%create probabibility distribution
for x=1:maxcodes
    Px(x) = zetafun(x,s);
end

%normalise Probability distribution
Px=Px./sum(Px); 

entropy=sum(-Px.*log2(Px));

codewords = createExpGcodes(maxcodes,k);
            
            codelengths = sum(~isnan(codewords),2);
            
          
            
            averagecodewordlength=sum(Px.*codelengths');
            
end