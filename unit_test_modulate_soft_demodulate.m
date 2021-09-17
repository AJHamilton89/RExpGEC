%% Unit test - Modulate (QPSK)

SNR = 100;
num_syms=10000;

TxCodeword= randi([0 1],num_syms,1)

%shape the bits
Txbits_QPSK = reshape(TxCodeword,[2,length(TxCodeword)/2]); %reshape the codeword to get ready to Tx - should be integrated into a function.

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

xhat = (sign(channelLLRs) + 1) / 2;
isequal(TxCodeword,xhat)

