% Copyright (C) 2021  Alex Hamilton

% This program is free software: you can redistribute it and/or modify it
% under the terms of the MIT Licence

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the MIT License for more details.



%%%%%%%%%%%%%%%%%
% inputs
%
% input bits - vector of bits of variable length (relating to the error
% corrected bits, post ExpGEC)
%

%
% outputs
%
% symbols - array of values corresponding to indexes in the input
% dictionary
%
%%%%%%%%%%%%%%%%%%

function [symbols] = ExpG_symbol_decoder(input,k,maxcodes)



% Determine number of codeword bits
% numbits = width(transitions)-4;

%initialise Look up table
maxwidth = length(de2bi(maxcodes+2^k-1))-1-k+length(de2bi(maxcodes+2^k-1));
codeset=NaN(maxcodes,maxwidth);


%this loop creates lookup table of all ExpG codes - LUT of codewords
for p=1:maxcodes
    codewordtemp=createExpGcode(p,k);
    codeset(p,1:length(codewordtemp))=codewordtemp;
end

codeset(isnan(codeset))=2; %bit of a bodge, but replace all the NaNs with a 2...



symbolnumber=1;%initialise symbol number
UECcount=0;%initalise as 0 - as the first bit will always be a UEC bit


n=1;

while n < length(input) %outer loop is searching for FLC states
    
    bit=input(n);
    
    if bit == 1 % move to FLC once we have found a 1 in the UEC codeword
        
        if n+UECcount+k > size(input,2) % error handling case for wrong size
            symbol = 0;
        else
            
            symbolbits = input(n-UECcount:n+UECcount+k);
            
            NANfiller = NaN(1,width(codeset)-length(symbolbits));
            row = [symbolbits NANfiller]; %make it look like a row in the codese
            row(isnan(row))=2;%make it look like a row in the codeset and do the same bodge as above
            
            if size(row,2) > maxwidth
                
                
                symbol = 0;
            else
                
                symbol = find(ismember(codeset, row,'rows'));%lookup in the codeset
                
            end
        end
        
        if isempty(symbol)%error handling case
            
            symbol = 0;
            
        end
        
        
        symbols(symbolnumber) = symbol;
        symbolnumber=symbolnumber + 1;%index +1
        
        
        
        n=n+UECcount+k+1;
        
        UECcount = 0;
        
    else
        
        UECcount = UECcount + 1; %count how many UEC bits have gone
        
        
        
        n=n+1;
        
    end
    
end

