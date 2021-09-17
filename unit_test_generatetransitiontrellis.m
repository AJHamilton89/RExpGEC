
%%%% - test script for trellis - checks that the above transitions screate
%%%% the  same as the bottom ones, by empirically testing many codewords

clear all


% set parameters

k=1;
depth=1; 
codingrate=2;
maxcodes=32;
s=3;
num_symbols=100000;


%Generate the trellis
trellis=generatetransitionstrellis(k,depth,codingrate); 

%%% Transmit Functionality
    
%Generate array of symbols in the zeta distribution
symbols=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);

%Generate a reordered ExpG codeword from the symbols
reorderedcodeword = generate_RExpGcodeword(k,symbols);

%work through trellis and return how many transitions were recorded
trellis=calculatetrellisprobs(trellis,reorderedcodeword);

transitions=trellis(:,(width(trellis)-1));

error= ( sum(transitions(1:(length(transitions)/2)))   -   sum(transitions((length(transitions)/2+1):(length(transitions))))  ) / sum (transitions)

