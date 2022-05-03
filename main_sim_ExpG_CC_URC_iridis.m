
%%%% - Test script for generating and transmitting RExpGEC codewords
function main_sim_ExpG_CC_URC_iridis(SNR_start, SNR_delta, SNR_stop, BER_stop, processes,p1,k,num_symbols,codingrate)





% set parameters


maxcodes = 1000; % set what the maximum value of codeset is. - Powers of 2  make sense

max_its=100;
num_its=max_its;
num_test_symbols=10000;
maxerrors = 10000; % this will affect how smooth the BER plots are 10000





%%%%% HPC functionality

% Set default values if input parameters are not provided
if ~exist('SNR_start','var')
    SNR_start = 10;
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
    p1 = 0.9;
end
if ~exist('k','var')
    k = 1;
end
if ~exist('depth','var')
   	depth = 1;
end
if ~exist('num_symbols','var')
    num_symbols = 100;
end
if ~exist('codingrate','var')
    codingrate = 2;
end

%Generate the trellis
% trellis=generatetransitionstrellis(k,depth,codingrate);

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


%% Complexity Calculation Section


% URC - this is the same as before

no_alphas = 2;
no_betas = 2;
no_gammas = 4;
no_deltas = 4;
av_codeword_length = Find_av_codeword_length_zeta(maxcodes,k,s);
no_trellis = av_codeword_length * num_symbols * codingrate;

Complexity_per_iter_URC = no_trellis*no_gammas+2*no_trellis*no_alphas+2*no_trellis*no_betas+2*no_trellis*no_deltas+3*no_trellis;


complexity_array_initial = 1:1:num_its;
complexity_array =  complexity_array_initial*Complexity_per_iter_URC;
%%


results = zeros(1,num_its+5);


% Choose a file to save the results into.
filename = ['VariablesStorage/ExpGCCURCresults_k=',num2str(k),'_R=',num2str(codingrate),'_p1=',num2str(p1),'_num_sym=',num2str(num_symbols),'_maxcode=',num2str(maxcodes),'_depth=',num2str(depth),'_snr=',num2str(SNR_start),'.mat'];
save(filename, 'results', '-MAT');

% Setup the SNR for the first iteration of the loop.
SNR_count = 1;
SNR = SNR_start;
FER = 1;
BER = 1;
BER1iter = 1;
BER10iter = 1;
BERRaw=1;

%% this section is for calculating the probabilities of the trellis - in practice this would be from an equation










%Generate array of symbols in the zeta distribution
symbolarray=generate_zeta_symbols_finite_dict(num_symbols*100,maxcodes,s);

codewordarray=generate_ExpGcodeword(k,symbolarray);

LLRcodewordones=log(sum(codewordarray)/length(codewordarray));

LLRcodewordzeroes=log([length(codewordarray)-sum(codewordarray)]/length(codewordarray));

LLRinputCC = LLRcodewordzeroes - LLRcodewordones ;


while SNR <= SNR_stop && BERRaw >= BER_stop
    
    % Counters to store the number of errors and bits simulated so far.
    errorbits=zeros(1,max_its);
    totalbits=0;
 
    rawerrors=0;
    totalrawbits=0;
    correctframes=0;
    totalframes=0;
        symbolsinerr=0;
    totalsymbols=0;
    
    while errorbits(num_its) < maxerrors && totalbits < (maxerrors/BER_stop) %if second logical statement then look at lower iteration count
        
        totalframes=totalframes+1;
        
        %%% Transmit Functionality
        startpoint=randi(num_symbols*99);
        symbols=symbolarray(startpoint:(num_symbols+startpoint));
        
        %Generate a reordered ExpG codeword from the symbols
        codeword = generate_ExpGcodeword(k,symbols);
        
        %Generate ExpG_CC codeword
        ExpG_CC=CC2_encoder(codeword);
        ExpG_CC=reshape(ExpG_CC,1,[]);
        
        %Generate the first interleaver & interleave
        interleaver1 = randperm(length(ExpG_CC)); % create the interleaver
        ExpG_CCint=ExpG_CC(interleaver1); % interleave the codeword
        
        %Pass through URC encorder
        TxCodeword=URC_encoder(ExpG_CCint);
        
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
        
        ExpG_CCinttildea = zeros(size(channelLLRs)); % vector = log(sum(codeword)/length(codeword))
        
        %ExpG_CCtildepa=zeros(1,size(channelLLRs,2)/2);
        
        ExpG_CCtildepa = repmat(LLRinputCC,1,(size(channelLLRs,2)/2));
        
        %% iterations will start from here
        for m=1:max_its
            %URC decoding
            
            ExpG_CCinttildee = URC2_decoder_bcjr(ExpG_CCinttildea,deinterleavedchannelLLRs); %
            
            
            
            % Deinterleaver 1
            ExpG_CCtildea = zeros(size(ExpG_CCinttildee));
            ExpG_CCtildea(interleaver1) = ExpG_CCinttildee;
            
             ExpG_CCtildea=reshape(ExpG_CCtildea,2,[]);
            
            % Trellis Decoder
            [ExpG_CCtildep,ExpG_CCtildee] = CC2_decoder_bcjr(ExpG_CCtildepa,ExpG_CCtildea);  %ExpG_CCtildepa should have a bias based on the probability of 1 vs 0 from source LLR of ratio.
            
            
            
            % Check if decoding has been successful
            xhat = (sign(ExpG_CCtildep) + 1) / 2;
            
            errorbits(m)=sum(codeword~=xhat) + errorbits(m);
            
            
            if isequal(xhat,codeword)
                correctframes=correctframes+1;
                break;
            end
            
            % Interleaver 1
            ExpG_CCinttildea = ExpG_CCtildee(interleaver1);
            ExpG_CCtildepa=ExpG_CCtildep;
        end
        
        totalbits=length(codeword) + totalbits;
        
        
        %%%%
%         Need to create a ExpG codeword decoder from CC decoder or
%         alternatively could do a packet/frame error rate?
%         
        %%%%
        
        
         Rxsymbols= ExpG_symbol_decoder(xhat,k,maxcodes);
         symbolsinerr = levenshtein_distance(Rxsymbols,symbols) + symbolsinerr;
         totalsymbols = length(symbols) + totalsymbols ;
        
        BER=errorbits(num_its)./totalbits;
%         BER10iter=errorbits(10)./totalbits;
%         BER1iter=errorbits(1)./totalbits;
        BERarray=errorbits/totalbits;
        SER=symbolsinerr/totalsymbols;
        FER=(totalframes-correctframes)/totalframes;
        BERRaw=rawerrors/totalrawbits;
        
        % Store the SNR and BER in a matrix and display it.
        results(SNR_count,1) = SNR;
        results(SNR_count,2) = errorbits(max_its);
        results(SNR_count,3) = totalbits;
        results(SNR_count,4) = FER;
        results(SNR_count,5) = SER;
        results(SNR_count,6:max_its+5) = BERarray;%BER10iter;
%         results(SNR_count,7) = BER1iter;
        
       
    end
     % Save the results into binary files. This avoids the loss of precision that is associated with ASCII files.
        save(filename, 'results', 'complexity_array' , '-MAT');
        
    
    if totalbits > (maxerrors/BER_stop)% this if statement will mean that from now on (going up in SNR) only the number of iterations that produce errors will be used.
        
        num_its=max(find(errorbits>maxerrors));
        
    end
    % Setup the SNR for the next iteration of the loop.
    SNR = SNR + SNR_delta;
    SNR_count = SNR_count + 1;
    
end
end







