
%%%% - Test script for generating and transmitting RExpGEC codewords
tic

clear all

% Channel SNR
snr = -5:0.01:5;

% set parameters
t1=clock;
plot_exit=0; %plot EXIT trajectories
k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 10000; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=100;
p1 = 0.555;
codingrate=2;
num_its=50;
num_test_symbols=10000;
num_runs=1;
s=zeta_p1_to_s(p1);
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

errorbits = zeros( length(snr),num_its);
totalbits = zeros (1, length(snr));


symbolsinerr = zeros (1, length(snr));
totalsymbols = zeros (1, length(snr));
BER = NaN(length(snr),num_its);
SER = NaN(1, length(snr));

figure (1)
xlabel('SNR [dB]');
ylabel('SER/BER')
ylim ([0.01 1])
grid minor

%Generate array of symbols in the zeta distribution
symbolarray=generate_zeta_symbols_finite_dict(num_symbols*100,maxcodes,s);


%%iridis HPC functionality

array = str2num(getenv('SLURM_ARRAY_TASK_ID'));


for index = 1:length(snr)
    
    
    while errorbits(index,(num_its)) < 10000
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
        N0 = 1/(10^(snr(index)/10));
        noise = sqrt(N0/2)*((randn(size(TxSignal)))+i*randn(size(TxSignal)));
        
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
            [RExpGECtildee, RExpGtildep] = RExpGEC_trellis_decoder(RExpGECtildea,trellis,probs);
            
            if plot_exit == 1
                %Measure Mutual Information
                MI_e_URC(m)=measure_mutual_information_averaging(RExpGECinttildee);
                MI_a_URC(m)=measure_mutual_information_averaging(RExpGECinttildea);
                MI_e_RExpGEC(m)=measure_mutual_information_averaging(RExpGECtildee);
                MI_a_RExpGEC(m)=measure_mutual_information_averaging(RExpGECtildea);
            end
            
            % Check if decoding has been successful
            xhat = (sign(RExpGtildep) + 1) / 2;
            
            errorbits(index,m)=sum(reorderedcodeword~=xhat) + errorbits(index,m);
            
            
            if isequal(xhat,reorderedcodeword)
                
                break;
            end
            
            % Interleaver 1
            RExpGECinttildea = RExpGECtildee(interleaver1);
        end
        
        %%% EXIT Charts
        if plot_exit == 1
            
            figure
            hold on
            axis square;
            ylabel('I_A');
            xlabel('I_E');
            xlim([0,1]);
            ylim([0,1]);
            
            plot(MI_a_URC,MI_e_URC,'d')
            plot(MI_e_RExpGEC,MI_a_RExpGEC,'d')
            stairs(MI_a_URC,MI_e_URC,'*-.')
            drawnow
            %%% EXIT Chart END
            
        end
        
        totalbits(index)=length(reorderedcodeword) + totalbits(index);
        
        Rxsymbols= RExpGEC_symbol_decoder(xhat,trellis,k,maxcodes);
        symbolsinerr(index) = levenshtein_distance(Rxsymbols,symbols) + symbolsinerr(index);%length(symbols)-sum(Rxsymbols == symbols);
        totalsymbols(index) = length(symbols) + totalsymbols(index) ;
        
        BER(index,:)=errorbits(index,:)./totalbits(index);
        SER(index)=symbolsinerr(index)/totalsymbols(index);
        
        
    end
    
    figure(1)
    hold on
    semilogy(snr,BER);
    semilogy(snr,SER,'LineWidth',2.0);
    
    ylim ([0.001 1])
    plottitle=sprintf('BER/SER vs SNR  Depth=%i Rate=%i K=%i maxcodes=%i num symbols=%i p1=%i.fig',depth,codingrate,k,maxcodes,num_symbols,p1);
                title(plottitle);
    grid minor
    drawnow
    set(gca, 'YScale', 'log')
    time=toc;
    
    fn1 = sprintf('Figures/BERSERCurveDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_s=%i_num_its=%i_SNR_start=%i_SNR_stop=%i.fig',depth,codingrate,k,maxcodes,num_symbols,s,num_its,snr(1),snr(index));
saveas(gcf,fn1)

fn2 = sprintf('VariablesStorage/BERSERWorkspaceDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_s=%i_num_its=%i.fig_SNR_start=%i_SNR_stop=%i',depth,codingrate,k,maxcodes,num_symbols,s,num_its,snr(1),snr(index));
save(fn2)
end








