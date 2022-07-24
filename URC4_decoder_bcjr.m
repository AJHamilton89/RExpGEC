% BCJR algorithm for an accumulator.
% Copyright (C) 2013  Robert G. Maunder

% apriori_uncoded_llrs is a 1xN vector of a priori uncoded LLRs
% apriori_encoded_llrs is a 1xN vector of a priori encoded LLRs
% extrinsic_uncoded_llrs is a 1xN vector of extrinsic encoded LLRs
function extrinsic_uncoded = URC4_decoder_bcjr(apriori_uncoded, apriori_encoded)

    if(length(apriori_uncoded) ~= length(apriori_encoded))
        error('LLR sequences must have the same length');
    end

    % Matrix to describe the trellis
    % Each row describes one transition in the trellis
    % Each state is allocated an index 1,2,3,... Note that this list starts
    % from 1 rather than 0.
    %               FromState,  ToState,    UncodedBit, EncodedBit
transitions =  [1,          1,          0,          0;
                1,          3,          1,          1;
                2,          3,          0,          1;
                2,          1,          1,          0;
                3,          4,          0,          1;
                3,          2,          1,          0;
                4,          2,          0,          0;
                4,          4,          1,          1];
               
    % Find the largest state index in the transitions matrix           
    % In this example, we have eight states since the code has three memory elements
    state_count = max(max(transitions(:,1)),max(transitions(:,2)));

    % Calculate the a priori transition log-probabilities
    gammas_encoded = zeros(size(transitions,1),length(apriori_uncoded));
    gammas_encoded(transitions(:,4)==1,apriori_encoded~=inf) = repmat(apriori_encoded(apriori_encoded~=inf), sum(transitions(:,4)==1),1);
    gammas_encoded(transitions(:,4)==0,apriori_encoded==inf) = -inf;
    
    % Calculate the a priori transition log-probabilities
    gammas = zeros(size(transitions,1),length(apriori_uncoded));
    gammas(transitions(:,3)==0 & transitions(:,4)==1,apriori_encoded~=inf) = repmat(apriori_encoded(apriori_encoded~=inf), sum(transitions(:,3)==0 & transitions(:,4)==1),1);
    gammas(transitions(:,3)==1 & transitions(:,4)==0,apriori_uncoded~=inf) = repmat(apriori_uncoded(apriori_uncoded~=inf), sum(transitions(:,3)==1 & transitions(:,4)==0),1);
    gammas(transitions(:,3)==1 & transitions(:,4)==1,apriori_uncoded~=inf&apriori_encoded~=inf) = repmat(apriori_uncoded(apriori_uncoded~=inf&apriori_encoded~=inf) + apriori_encoded(apriori_uncoded~=inf&apriori_encoded~=inf), sum(transitions(:,3)==1 & transitions(:,4)==1),1); % 1 addition
    gammas(transitions(:,3)==0 & transitions(:,4)==0,apriori_uncoded==inf|apriori_encoded==inf) = -inf;

    % Recursion to calculate forward state log-probabilities
    alphas=zeros(state_count,length(apriori_uncoded));
    alphas(2:end,1)=-inf; % We know that these are not the first state
    for bit_index = 2:length(apriori_uncoded)        
        temp = alphas(transitions(:,1),bit_index-1)+gammas(:,bit_index-1);
        for state_index = 1:state_count
            alphas(state_index,bit_index) = maxstar(temp(transitions(:,2) == state_index));
        end
    end
    
    % Recursion to calculate backward state log-probabilities
    betas=zeros(state_count,length(apriori_uncoded));
    for bit_index = length(apriori_uncoded)-1:-1:1
        temp = betas(transitions(:,2),bit_index+1)+gammas(:,bit_index+1);
        for state_index = 1:state_count
            betas(state_index,bit_index) = maxstar(temp(transitions(:,1) == state_index));
        end
    end

    % Calculate a posteriori transition log-probabilities
    deltas = alphas(transitions(:,1),:) + betas(transitions(:,2),:) + gammas_encoded;
    
    % Calculate the uncoded extrinsic LLRs
    log_p0=maxstar(deltas(transitions(:,3) == 0,:));
    log_p1=maxstar(deltas(transitions(:,3) == 1,:));
    extrinsic_uncoded = log_p1-log_p0;
   
