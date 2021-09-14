
%%%% - Test script for generating and transmitting RExpGEC codewords



% set parameters

k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 32; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=10;
s = 3;
codingrate=2;
numruns = 100;
SNR=0;
num_its=100;
num_test_symbols=10000;

% calculate dependent parameters


%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate);

%noise power spectral density
N0 = 1/(10^(SNR/10));

%% this section is for calculating the probabilities of the trellis - in practice this would be from an equation

%Generate array of symbols in the zeta distribution
symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);

%Generate a reordered ExpG codeword from the symbols
reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);


%work through trellis and return how many transitions were recorded
probs=calculatetrellisprobs(trellis,reorderedcodeword_probs);

%Array to check if it worked
Positive_result=zeros(2,numruns);




for n=1:numruns
    
    %%% Transmit Functionality
    
    %Generate array of symbols in the zeta distribution
    symbols=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);
    
    %Generate a reordered ExpG codeword from the symbols
    reorderedcodeword = generate_RExpGcodeword(k,symbols);
    
    %Generate RExpGEC codeword
    RExpGEC=generateRExpGEC(trellis,reorderedcodeword);
    
    %Generate the first interleaver & interleave
    interleaver1 = randperm(length(RExpGEC)); % create the interleaver
    RExpGECint=RExpGEC(interleaver1); % interleave the codeword
    
    %Pass through URC encorder
    TxCodeword=URC_encoder(RExpGECint);
    
    %Generate the second interleaver & interleave
    interleaver2 = randperm(length(TxCodeword)); % create the interleaver
    TxCodewordint=TxCodeword(interleaver2); % interleave the codeword
    
    %shape the bits
    Txbits_QPSK = reshape(TxCodewordint,[2,length(TxCodewordint)/2]); %reshape the codeword to get ready to Tx - should be integrated into a function.
    
    %Modulate the codeword
    TxSignal = modulate(Txbits_QPSK);
    
    
    
    %create the channel
    channel = sqrt(1/2) * (randn(size(TxSignal)) + 1i*randn(size(TxSignal)));
    
    %Create the Rx Signal simulated from channel
    RxSignal = channel.*TxSignal + sqrt(N0/2) * (randn(size(TxSignal)) + 1i*randn(size(TxSignal)));
    
    
    
    %%% Receive Functionality
    
    %demodulate
    LLRs_QPSK_tilde = soft_demodulate(RxSignal, channel, N0);
    channelLLRs = reshape(LLRs_QPSK_tilde,size(TxCodeword));
    
    %deinteleave - interleaver 2
    deinterleavedchannelLLRs = zeros(size(channelLLRs));
    deinterleavedchannelLLRs(interleaver2) = channelLLRs;
    
    RExpGECinttildea = zeros(size(channelLLRs));
    
    %% iterations will start from here
    for m=1:num_its
        %URC decoding
        
        RExpGECinttildee = URC2_decoder_bcjr(RExpGECinttildea,deinterleavedchannelLLRs); %
        
        % Deinterleaver 1
        RExpGECtildea = zeros(size(RExpGECinttildee));
        RExpGECtildea(interleaver1) = RExpGECinttildee;
        
        % Trellis Decoder
        [RExpGECtildee, RExpGtildep] = RExpGEC_trellis_decoder(RExpGECtildea,trellis,probs,num_symbols);
        
        % Check if decoding has been successful
        xhat = (sign(RExpGtildep) + 1) / 2;;
        if isequal(xhat,reorderedcodeword)
            Positive_result(1,n)=1;
            Positive_result(2,n)=m;
            break;
        end
        
        % Interleaver 1
        RExpGECinttildea = RExpGECtildee(interleaver1);
    end
    
end

