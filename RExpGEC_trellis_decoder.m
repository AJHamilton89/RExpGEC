% Copyright (C) 2013  Robert G. Maunder

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
function [RExpGECtildee, RExpGtildep] = RExpGEC_trellis_decoder(RExpGECtildea,trellis,probs,a)


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


% Determine number of transitions
r = size(trellis,1);


% Determine if ending at same side of trellis or opposite
if mod(a,2) ==0
    endstateaddition = 0;
else
    endstateaddition = r/2;
end

% Determine number of transitions
maxstate = max(trellis(:,1));

% Determine number of codeword bits
n = width(trellis)-4;

% Reshape codeword - is this needed - reshapes into an array the depth of
% each EC code
RExpGECtildea=reshape(RExpGECtildea,n,[]);

% Determine length of unary-encoded bit sequence
b = size(RExpGECtildea,2);

% Build the trellis and calculate conditional transition probabilities
transitions = trellis;
transitions(:,4) = [] ;%remove the column relating to  leaf nodes as it is superfluous

%create non-zero states
transitions(:,1:2)= transitions(:,1:2)+1;

%set gammas up from the probs
gammas = repmat(log(probs),1,b);

% Calculate a priori transition log-probabilities
for codebit_index = 1:n
    gammas(transitions(:,3+codebit_index)==1,:) = gammas(transitions(:,3+codebit_index)==1,:) + repmat(RExpGECtildea(codebit_index, :), sum(transitions(:,3+codebit_index)==1),1);
end

% Forward recursion to calculate state log-probabilities
alphas=-inf(r,b);
alphas(1,1)=0; % We know that this is the first state
for j = 2:b
    temp = alphas(transitions(:,1),j-1)+gammas(:,j-1);
    for m = 1:maxstate
        alphas(m,j) = maxstar(temp(transitions(:,2) == m));
    end
end

% Backward recursion to calculate state log-probabilities
betas=-inf(r,b);
betas(1 + endstateaddition,end)=0; % We know that this is the final state - (or the other one...)we don't for RExpG...
for j = b-1:-1:1
    temp = betas(transitions(:,2),j+1)+gammas(:,j+1);
    for mprime = 1:maxstate
        betas(mprime,j) = maxstar(temp(transitions(:,1) == mprime));
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
RExpGECtildee = yztildep(2:end,:) -RExpGECtildea;

%reshape in order to pass back
RExpGECtildee = reshape(RExpGECtildee,1,[]);

end