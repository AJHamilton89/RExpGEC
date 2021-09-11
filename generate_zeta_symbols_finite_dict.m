%inputs - num_symbols - how many symbols you want; s- s parameter of zeta
%function
% output is an array of symbol values

function symbols = generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s)


%create probabibility distribution
for x=1:maxcodes
    Px(x) = zetafun(x,s);
end

%normalise Probability distribution
Px=Px./sum(Px); 

% generate symbols
symbols=generate_random_symbols_Px(num_symbols,Px); %

end
