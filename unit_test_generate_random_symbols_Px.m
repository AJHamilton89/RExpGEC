% unit test for the function generate_random_symbols_Px
% Alex Hamilton - 11/9/2021

num_runs=10;

max_sym = 100000;
min_sym = 100;

min_dist=20;
max_dist=1000;

for n=1:num_runs

num_symbols_to_test= round(min_sym + (max_sym-min_sym) .* rand); 
size_of_distribution= round(min_dist + (max_dist-min_dist) .* rand);


%generate a probability distribution (random)
Px = rand(1,size_of_distribution);
Px=Px./sum(Px);


%generate the symbols
symbols=generate_random_symbols_Px(num_symbols_to_test,Px);

%calculate the histogram
histbins=histc(symbols,0.5:1:(size_of_distribution-0.5));
normalisedhist=histbins./sum(histbins);

%calculate the error accross all symbols
errorvector(n)=sum(Px-normalisedhist)/size_of_distribution;

standardeviation(n)=errorvector(n)/sqrt(size_of_distribution);

end
