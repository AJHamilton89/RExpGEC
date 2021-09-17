
%%%% - Calculating the trellis probabilities via simulation

clear all


% set parameters

k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 32; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=10000;
s = 3;
codingrate=2;
numruns = 1;

%% section only used for debugging
%create probabibility distribution
for x=1:maxcodes
    Px(x) = zetafun(x,s);
end

%normalise Probability distribution
Px=Px./sum(Px); 
%%

%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate); 



%%% Transmit Functionality
    
%Generate array of symbols in the zeta distribution
symbols=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);

%Generate a reordered ExpG codeword from the symbols
reorderedcodeword = generate_RExpGcodeword(k,symbols);


%work through trellis and return how many transitions were recorded
[probs,trellis]=calculatetrellisprobs(trellis,reorderedcodeword);
