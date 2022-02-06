
%%%% - Test script for generating and transmitting RExpGEC codewords
function main_sim_RExpGEC_iridis(SNR_start, SNR_delta, SNR_stop, BER_stop, processes,p1)





% set parameters


k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 1000; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=1000;

codingrate=2;
num_its=50;
num_test_symbols=10000;
maxerrors = 10000; % this will affect how smooth the BER plots are 10000


%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate);

%%%%% HPC functionality

% Set default values if input parameters are not provided
if ~exist('SNR_start','var')
    SNR_start = -10;
end
if ~exist('SNR_delta','var')
    SNR_delta = 0.25;
end
if ~exist('SNR_stop','var')
    SNR_stop = inf;
end
if ~exist('BER_stop','var')
    BER_stop = 0;
end
if ~exist('processes','var')
    processes = 1;
end
if ~exist('p1','var')
    p1 = 0.555;
end

% Deal with parallel processing on slurm
process = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(process)
    process = 0;
    processes = 1;
else
    if process >= processes
        error('process >= processes');
    end
end
SNR_start = SNR_start+SNR_delta*process;
SNR_delta = SNR_delta*processes;


s=zeta_p1_to_s(p1);
% calculate dependent parameters



results = zeros(1,num_its+5);

% Choose a file to save the results into.
filename = ['VariablesStorage/results_k=',num2str(k),'_R=',num2str(codingrate),'_p1=',num2str(p1),'_num_sym=',num2str(num_symbols),'_snr=',num2str(SNR_start),'.mat'];
save(filename, 'results', '-MAT');

% Setup the SNR for the first iteration of the loop.
SNR_count = 1;
SNR = SNR_start;
BER = 1;
BER1iter = 1;
BER10iter = 1;
BERRaw=1;

%% this section is for calculating the probabilities of the trellis - in practice this would be from an equation

%Generate array of symbols in the zeta distribution
symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);

%Generate a reordered ExpG codeword from the symbols
reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);


%work through trellis and return how many transitions were recorded
probs=calculatetrellisprobs(trellis,reorderedcodeword_probs);








%Generate array of symbols in the zeta distribution
symbolarray=generate_zeta_symbols_finite_dict(num_symbols*100,maxcodes,s);





while SNR <= SNR_stop && BERRaw >= BER_stop
    
    % Counters to store the number of errors and bits simulated so far.
    errorbits=zeros(1,num_its);
    totalbits=0;
    symbolsinerr=0;
    totalsymbols=0;
    rawerrors=0;
    totalrawbits=0;
    
    while errorbits(num_its) < maxerrors && totalbits < (maxerrors/BER_stop) %if second logical statement then look at lower iteration count
        
        
        
        %%% Transmit Functionality
        startpoint=randi(num_symbols*99);
        symbols=symbolarray(startpoint:(num_symbols+startpoint));
        
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
        %channel = sqrt(1/2)*(randn(size(TxSignal)))+i*randn(size(TxSignal));
        
        % AWGN channel
        channel = ones(size(TxSignal));
        
        % Generate some noise
        N0 = 1/(10^(SNR/10));
        noise = sqrt(N0/2)*((randn(size(TxSignal)))+i*randn(size(TxSignal)));
        
        % Generate the received signal
        RxSignal = TxSignal.*channel+noise;
        
        
        
        
        %%% Receive Functionality
        
        %demodulate
        LLRs_QPSK_tilde = soft_demodulate(RxSignal, channel, N0);
        channelLLRs = reshape(LLRs_QPSK_tilde,size(TxCodeword));
        
        rawxhat = (sign(channelLLRs) + 1) / 2;
        rawerrors=sum(TxCodewordint~=rawxhat) + rawerrors;
        
        totalrawbits=length(TxCodewordint) + totalrawbits;
        
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
            [RExpGECtildee, RExpGtildep] = RExpGEC_trellis_decoder(RExpGECtildea,trellis,probs);
            
            
            
            % Check if decoding has been successful
            xhat = (sign(RExpGtildep) + 1) / 2;
            
            errorbits(m)=sum(reorderedcodeword~=xhat) + errorbits(m);
            
            
            if isequal(xhat,reorderedcodeword)
                
                break;
            end
            
            % Interleaver 1
            RExpGECinttildea = RExpGECtildee(interleaver1);
        end
        
        totalbits=length(reorderedcodeword) + totalbits;
        
        Rxsymbols= RExpGEC_symbol_decoder(xhat,trellis,k,maxcodes);
        symbolsinerr = levenshtein_distance(Rxsymbols,symbols) + symbolsinerr;
        totalsymbols = length(symbols) + totalsymbols ;
        
        BER=errorbits(num_its)./totalbits;
        BER10iter=errorbits(10)./totalbits;
        BER1iter=errorbits(1)./totalbits;
        BERarray=errorbits/totalbits;
        SER=symbolsinerr/totalsymbols;
        BERRaw=rawerrors/totalrawbits;
        
        % Store the SNR and BER in a matrix and display it.
        results(SNR_count,1) = SNR;
        results(SNR_count,2) = errorbits(num_its);
        results(SNR_count,3) = totalbits;
        results(SNR_count,4) = BER;
        results(SNR_count,5) = SER;
        results(SNR_count,6:num_its+5) = BERarray;%BER10iter;
%         results(SNR_count,7) = BER1iter;
        
        % Save the results into binary files. This avoids the loss of precision that is associated with ASCII files.
        save(filename, 'results', '-MAT');
        
    end
    
    if totalbits > (maxerrors/BER_stop)% this if statement will mean that from now on (going up in SNR) only the number of iterations that produce errors will be used.
        
        num_its=max(find(errorbits>maxerrors));
        
    end
    % Setup the SNR for the next iteration of the loop.
    SNR = SNR + SNR_delta;
    SNR_count = SNR_count + 1;
    
end
end







