% function to generate RExpG codeword based upon k parameter, and an input
% array of the symbol values that wish to be transmitted.
% inputs k = k parameter; sybols = array of symbols


function codeword = generate_ExpGcodeword(k,symbols)

%initalise codeword
codeword=[];

%infer the maximum codes required
maxcodes=max(symbols);

%initialise Look up table
codeset=NaN(maxcodes,100);

%this loop creates lookup table of all RExpG codes - LUT of codewords
for p=1:maxcodes
    codewordindex=createExpGcode(p,k);
    codeset(p,1:length(codewordindex))=codewordindex;
end

%create reorded codeword string 
for n=1:length(symbols)
    codewordtemp = codeset (symbols(n),:);
    codewordtemp=(codewordtemp(~isnan(codewordtemp)));
    codeword = [codeword,codewordtemp];
    
end

end
