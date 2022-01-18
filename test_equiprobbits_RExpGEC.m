% Script for testing equiprobable bits
% Copyright (C) 2021  Alex Hamilton
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% The GNU General Public License can be seen at http://www.gnu.org/licenses/.
clear all

% set REXpGEC parameters

%variables that don't loop



maxcodes = 100; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=100; %defines block size
p1_count=10;
p1 = 0.00001+0.9999*(0:1/(p1_count-1):1);

%set some parameters

block_count = 1000;
maxbits=100;
varloop=1;
totalbits=0;
RExpGECstorage=zeros(p1_count,maxbits);
codingrate=2;
depth=1;
k=1;



%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate);


for p1_index = 1:p1_count
    
    s=zeta_p1_to_s(p1(p1_index));
    s_storage(p1_index)=s;
    
    %Generate array of symbols in the zeta distribution
    symbolsarray=generate_zeta_symbols_finite_dict(num_symbols*block_count,maxcodes,s);
    
    
    
    for block_index = 1:block_count
        
        symbols=symbolsarray(((block_index-1)*num_symbols+1):((block_index)*num_symbols));
        
        
        %Generate a reordered ExpG codeword from the symbols
        clear reorderedcodeword
        reorderedcodeword = generate_RExpGcodeword(k,symbols);
        
        %Generate RExpGEC codeword
        clear RExpGEC
        RExpGEC=generateRExpGEC(trellis,reorderedcodeword);
        
        lengthRExpGEC(varloop)=length(RExpGEC);
        
        sizemax=max(size(RExpGEC,2),size(RExpGECstorage(p1_index,:) ,2));
        sizemin=min(size(RExpGEC,2),size(RExpGECstorage(p1_index,:) ,2));
        
        RExpGECstoretemp = zeros (2,sizemax);
        
        RExpGECstoretemp(1,1:length(RExpGEC))=RExpGEC;
        RExpGECstoretemp(2,1:length(RExpGECstorage))=RExpGECstorage(p1_index,1:maxbits);
        
        RExpGECstorage (p1_index,1:maxbits) = RExpGECstoretemp(1,1:maxbits)+RExpGECstoretemp(2,1:maxbits);
        
        totalbits=sizemax+totalbits;
        
        varloop=varloop+1;
    end
    Bitprobs(p1_index,1:maxbits) = RExpGECstorage (p1_index,1:maxbits) / block_count;
end

for n=1:p1_count
    
    figure
    xlabel('bit location')
    ylabel('P(b=1)')
    bar(Bitprobs(n,:))
     plottitle=sprintf('Probability of bits in a location being 1, Rate=%i K=%i  num symbols=%i p1=%i.fig',codingrate,k,num_symbols,p1(p1_index));
   title(plottitle)
end


