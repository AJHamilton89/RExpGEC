% Copyright (C) 2021  Alex Hamilton

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at
% http://www.gnu.org/licenses/.

%%%%%%%%%%%%%%%%%
% inputs
%
% input bits - vector of bits of variable length (relating to the error
% corrected bits, post RExpGEC)
%
% trellis - defines the RExpG trellis
%
% outputs
%
% symbols - array of values corresponding to indexes in the input
% dictionary
%
%%%%%%%%%%%%%%%%%%

function [symbols] = RExpG_symbol_decoder(input,trellis,k,maxcodes)



% Determine number of codeword bits
% numbits = width(transitions)-4;

%initialise Look up table
maxwidth = length(de2bi(maxcodes+2^k-1))-1-k+length(de2bi(maxcodes+2^k-1));
codeset=NaN(maxcodes,maxwidth);


%this loop creates lookup table of all RExpG codes - LUT of codewords
for p=1:maxcodes
    codewordtemp=createExpGcode(p,k);
    reordererdtemp=reordercodewordRExpG(codewordtemp,k);
    codeset(p,1:length(reordererdtemp))=reordererdtemp;
end

codeset(isnan(codeset))=2; %bit of a bodge, but replace all the NaNs with a 2...

codewordtemp=[]; %initialise RExpGEC codeword as R
startstate=0;
fromstate=startstate;
symbolnumber=1;%initialise symbol number

for n=1:length(input)
    
    bit=input(n);
    
    possiblestatesid = trellis(:,3) == bit; %create subset of transitions matching the bit
    possiblestates = trellis(possiblestatesid,:);
    
    actualtransitionid = possiblestates(:,1) == fromstate; %navigate actual state
    actualtransition = possiblestates(actualtransitionid,:);
    
    tostate=actualtransition(2);
    
    codewordtemp = [codewordtemp bit];
    
    if (0<=tostate) && (tostate<=1) % if the trellis goes back to one of the start states
        
        NANfiller = NaN(1,width(codeset)-length(codewordtemp));
        row = [codewordtemp NANfiller]; %make it look like a row in the codese
        row(isnan(row))=2;%make it look like a row in the codeset and do the same bodge as above
        
        if length(row) > maxwidth %error handler for when the received 'symbol' is too long.
            
            symbols(symbolnumber) = 0;
            
        else
        
        symbol = find(ismember(codeset, row,'rows'));%lookup in the codeset
        
        
        
        
        
        if isempty(symbol)%error handling case - if we end up at stop point of trellis but symbol isn't found allocate symbol value of 0 so that it will be registered as an error
            
            symbols(symbolnumber) = 0;
            
        else
            
            symbols(symbolnumber) = symbol;
            
        end
        end
        symbolnumber=symbolnumber+1;%index +1
        
        codewordtemp= []; %reset codewordtemp
    end
    
    
    
    
    %     codeword(n*(width(trellis)-4)-(width(trellis)-5):n*(width(trellis)-4))=actualtransition(5:end); %populate RExpGEC - this should be paramaterisable
    
    
    fromstate=tostate; %create new from state to go back through loop
    
    
    
    
    
    
end

end

