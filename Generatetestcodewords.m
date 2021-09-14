 clear all
depth = 2;
maxcodes = 32; % set what the maximum value of codeset is. - Powers of 2  make sense
k=0;
codewordout=[];


for num_symbols=logspace(5,10,6);
s = 3;

for x=1:maxcodes
Px(x) = zetafun(x,s);
end

%normalise Probability distribution
Px=Px./sum(Px);

% 
% codewords = createExpGcodes(maxcodes,k);
treematrix = generatetransitions(depth,k);
%reorderedcodeword = zeros(size(codewords));

% for n=1:maxcodes
% codeword=idcodeword(n,codewords);
% reorderedcodeword=reordercodewordRExpG(codeword,k);
% codewordout = [codewordout,reorderedcodeword]; %currently produces a long string of codewords
% end

% generate symbols
symbols=generate_random_symbols(num_symbols,'distribution_zeta',s,maxcodes);


codeset=NaN(maxcodes,100);
for p=1:maxcodes;
codewordtemp=createExpGcode(p,k);
reordererdtemp=reordercodewordRExpG(codewordtemp,k);
codeset(p,1:length(reordererdtemp))=reordererdtemp;
end

for n=1:length(symbols)
% clear reorderedcodeword
 reorderedcodeword = codeset (symbols(n),:);
 reorderedcodeword=(reorderedcodeword(~isnan(reorderedcodeword)));

 % codeword=createExpGcode(symbols(n),k);
% 
% reorderedcodeword=reordercodewordRExpG(codeword,k);

codewordout = [codewordout,reorderedcodeword];

end

Filename = sprintf('codeword_%s_K=%d_S=%d_size=%d.mat', datestr(now,'mm-dd-yyyy HH-MM'),k, s, num_symbols);
save(Filename,'codewordout')
end
