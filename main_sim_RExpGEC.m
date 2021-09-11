
%%%% - Test script for generating and transmitting RExpGEC codewords



% set parameters

k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 32; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=10;
s = 3;
codingrate=2;
numruns = 1;

%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate); 

for n=1:numruns

%%% Transmit Functionality
    
%Generate array of symbols in the zeta distribution
symbols=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);

%Generate a reordered ExpG codeword from the symbols
reorderedcodeword = generate_RExpGcodeword(k,symbols);

%Generate RExpGEC codeword
RExpGEC=generateRExpGEC(trellis,reorderedcodeword);

%Modulate the codeword
TxSignal=QPSKmodulate(RExpGEC);% unit test against demodulate

%%% Receive Functionality

end

