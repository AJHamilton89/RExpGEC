%%%code to calcuate probabilties in accordance with page 59 of logbook

createcodes=1;
for max=[10];%10,100,1000,

for k=0:6
    index=1;
   
    
    
    if createcodes==1
codewords = createExpGcodes(max,k);
end
for p1=0.025:0.05:1
  
%set how many symbols we need to go down (also includes the test length)
num_symbols = 1000;

%set S parameter - only required for zeta
s=zeta_p1_to_s(p1);

%set K parameter for codewords


Pac=0;
Pab=0;
Pcf=0;
Pce=0;
Pdg=0;
Pgi=0;
Pgh=0;
Pdb=0;
Pbdzeros=0;
Pbdones=0;

%set flag for creating codewords


createsymbols=1;


%set how far we want to go for conditional probabilities
g=20;


%% create RExp codewords


%% select codeword if required
%identify symbol of interest

if createsymbols==1
%generate symbols
symbols=generate_random_symbols(num_symbols,'distribution_zeta',s,max);
end

codewordout=[];

for num=1:length(symbols)
codeword=idcodeword(symbols(num),codewords);
codewordout = [codewordout,codeword];
end
averagecodewordlength=length(codewordout)/length(symbols);



%geometric
% for x=1:max
% Px (x)= p1*(1-p1).^(x-1);
% end
%zeta

for x=1:max
Px (x) = x.^-s/zeta(s);
end

information=measure_information_generic(Px);

ratenew(index,k+1)=information/averagecodewordlength;
index=index+1;
end
end
figure


plot(0.025:0.05:1,ratenew)

ylabel('\eta_0')
title(['Information Efficiency vs P1 of zeta distributions L =' num2str(max)])
xlabel('P1')
ylim([0 1])
xlim([0 1])
grid on
drawnow
end