% Copyright (C) 2013  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

% The soft-input hard-output RExpG 
% ytildep is the LLR vector comprising b a posteriori LLRs, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
% a is the number of symbols to decode
% xhat is the symbol vector comprising a symbols
function xhat = RExpG_decoder_hard(RExpGtildep, a, maxcodes,k)



%initialise Look up table
codeset=NaN(maxcodes,100);

% create LUT of codewords
for p=1:maxcodes
    codewordtemp=createExpGcode(p,k);
    reordererdtemp=reordercodewordRExpG(codewordtemp,k);
    codeset(p,1:length(reordererdtemp))=reordererdtemp;
end

yhat = (sign(RExpGtildep) + 1) / 2; % take a hard decision

for n=1:size(yhat,2)
    
end
[~,indices] = sort(RExpGtildep,'descend');
yhat(indices(1:a)) = 1;
xhat = unary_decoder_hard(yhat);