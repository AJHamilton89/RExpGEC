
%%%% - Test script for generating and transmitting RExpGEC codewords

clear all

% Channel SNR
snr = 100:1:100;

% set parameters

k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 32; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=2;
s = 3;
codingrate=2;
num_its=10;
num_test_symbols=10000;
num_runs=10;

% calculate dependent parameters


%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate);



%% this section is for calculating the probabilities of the trellis - in practice this would be from an equation

%Generate array of symbols in the zeta distribution
symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);

%Generate a reordered ExpG codeword from the symbols
reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);


%work through trellis and return how many transitions were recorded
probs=calculatetrellisprobs(trellis,reorderedcodeword_probs);


Positive_result = zeros(1,num_runs);


for index = 1:length(snr)
    
    
    for n=1:num_runs
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
    
    
    
    
    
    % Channel
    % -------
    
    % Uncorrelated Rayleigh fading channel
    channel = sqrt(1/2)*(randn(size(TxSignal)))+i*randn(size(TxSignal));
    
    % AWGN channel
    %channel = ones(size(TxSignal));
   
    % Generate some noise
    N0 = 1/(10^(snr(index)/10));
    noise = sqrt(N0/2)*(randn(size(TxSignal)))+i*randn(size(TxSignal));
    
    % Generate the received signal
    RxSignal = TxSignal.*channel+noise;
    

    
    
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
            Positive_result(n)=1;
            break;
        end
        
        % Interleaver 1
        RExpGECinttildea = RExpGECtildee(interleaver1);
    end
    
    errorbits(index,n)=sum(reorderedcodeword~=xhat);
    totalbits(index,n)=length(reorderedcodeword);
    
  
    
    end
    
      Positive_result
end

errorbits2=sum(errorbits,2);
totalbits2=sum(totalbits,2);

BER=errorbits2./totalbits2;

%figure

% semilogy(snr,BER);
% xlabel('SNR [dB]');
% ylabel('BER')