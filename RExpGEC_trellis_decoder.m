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



% The soft-input soft-output trellis decoder of Section IV-A in (Maunder et al., 2013)
% ztildea is the LLR matrix comprising n times b a priori LLRs pertaining to z, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
% Px_vector is a vector comprising the probabilities of the first r/2-1 symbol values
% a is the number of symbols to decode
% ztildee is the LLR matrix comprising n times b extrinsic LLRs pertaining to z, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
% ytildep is the LLR vector comprising b a posteroiri LLRs pertaining to y, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
function [RExpGECtildee, RExpGtildep] = RExpGEC_trellis_decoder(RExpGECtildea,transitions,probs)


% All calculations are performed in the logarithmic domain in order to
% avoid numerical issues. These occur in the normal domain, because some of
% the confidences can get smaller than the smallest number the computer can
% store.
%
% A multiplication of two confidences is achieved using the addition of the
% corresponding log-confidences. If A = log(a) and B = log(b), then
% log(a*b) = A+B.
%
% An addition of two confidences is achieved using the Jacobian logarithm
% of the corresponding log-confidences. The Jacobian logarithm is defined
% in the jac.m file. If A = log(a) and B = log(b), then
% log(a+b) = max(A,B) + log(1+exp(-abs(A-B))).

% Determine number of codeword bits
n = width(transitions)-4;


% Reshape codeword - reshapes into an array the depth of
% each RExpG symbol / RExpGEC codeword
RExpGECtildea=reshape(RExpGECtildea,n,[]);

%create non-zero states on the transitions - because the zero indexed states don't work
transitions(:,1:2)= transitions(:,1:2)+1;

transitions(:,4) = [];

% Determine number of states
r = max(transitions(:,2));

% Determine length of RExpG encoded bit sequence
b = size(RExpGECtildea,2);

mid = round(b/2);

%set gammas up from the probs

if nargin < 3 %if we don't have the probs passed through
    
    gammas = zeros (size (transitions,1),b);
    
else

gammas = repmat(log(probs),1,b); % 
% gammas = repmat(log((probs./(1-probs))),1,b);

end

% Calculate a priori transition log-probabilities


for codebit_index = 1:n
    gammas(transitions(:,3+codebit_index)==1,:) = gammas(transitions(:,3+codebit_index)==1,:) + repmat(RExpGECtildea(codebit_index, :), sum(transitions(:,3+codebit_index)==1),1); 
%     gammas(transitions(:,3+codebit_index)==0,:) = gammas(transitions(:,3+codebit_index)==0,:) - repmat(RExpGECtildea(codebit_index, :), sum(transitions(:,3+codebit_index)==0),1);
end

%             complexity_counter = complexity_counter+size(transitions,1);


% Forward recursion to calculate state log-probabilities
alphas=-inf(r,b);
alphas(1,1)=0; % We know that this is the first state 


for j = 2:b
    temp = alphas(transitions(:,1),j-1)+gammas(:,j-1);                       
    for m = 1:r
        if ismember(m,transitions(:,2))
        alphas(m,j) = maxstar(temp(transitions(:,2) == m));
        end
    end
end
% 
%             complexity_counter = 2*sum(~isinf(alphas(transitions(:,1),mid)+gammas(:,mid))) + complexity_counter;

% Backward recursion to calculate state log-probabilities
betas=-inf(r,b);
betas(1:2,end)=0; % one of the end states 
    
for j = b-1:-1:1
    temp = betas(transitions(:,2),j+1)+gammas(:,j+1);
    for mprime = 1:r
        if ismember(mprime,transitions(:,1))
        betas(mprime,j) = maxstar(temp(transitions(:,1) == mprime));
        end
    end
end

% Calculate a posteriori transition log-probabilities
deltas = alphas(transitions(:,1),:) + betas(transitions(:,2),:) + gammas;

% Calculate a posteriori LLRs pertaining to y and z
yztildep = zeros(n+1,b);
for codebit_index = 1:n+1
    log_p1=maxstar(deltas(transitions(:,2+codebit_index) == 1,:));
    log_p0=maxstar(deltas(transitions(:,2+codebit_index) == 0,:));
    yztildep(codebit_index,:) = log_p1-log_p0;
end



% Extract a posteriori LLRs pertaining to y
RExpGtildep = yztildep(1,:);

% Calculate extrinsic LLRs pertaining to z
RExpGECtildee = yztildep(2:end,:) - RExpGECtildea;

%reshape in order to pass back
RExpGECtildee = reshape(RExpGECtildee,1,[]);

end