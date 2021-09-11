% function to generate RExpG codeword based upon k parameter, and an input
% array of the symbol values that wish to be transmitted.
% inputs k = k parameter; sybols = array of symbols


function reorderedcodeword = generate_RExpGcodeword(k,symbols)

%initalise codeword
reorderedcodeword=[];

%infer the maximum codes required
maxcodes=max(symbols);

%initialise Look up table
codeset=NaN(maxcodes,100);

%this loop creates lookup table of all RExpG codes - LUT of codewords
for p=1:maxcodes
    codewordtemp=createExpGcode(p,k);
    reordererdtemp=reordercodewordRExpG(codewordtemp,k);
    codeset(p,1:length(reordererdtemp))=reordererdtemp;
end

%create reorded codeword string 
for n=1:length(symbols)
    reorderedcodewordtemp = codeset (symbols(n),:);
    reorderedcodewordtemp=(reorderedcodewordtemp(~isnan(reorderedcodewordtemp)));
    reorderedcodeword = [reorderedcodeword,reorderedcodewordtemp];
    
end

end
