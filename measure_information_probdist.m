function information = measure_information_probdist(Px)




for n=1:length(Px)
    entropy(n)=-Px(n)*log2(Px(n));
end




information=sum(entropy(~isnan(entropy)));