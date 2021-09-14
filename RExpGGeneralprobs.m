%Generalise Test probabilities

k=2;
s = 2;
totalprobs=0;
totaltransitions=0;
start=3;
for n=1:k*2
    startdenomtemp(n)=zetafun(n,s);
end
startdenom=sum(startdenomtemp);
startnumerator=startdenom;
stop = 10000;
onetrans = 0;
zerotrans = 0;

 zerotranstemp=zeros(stop,1);
 onetranstemp=zeros(stop,1);
prob=zeros(stop,1);
normprob=zeros(stop,1);
 bival = createExpGcodes(stop,k);
 
 PAC = ( zetafun((2*k),s) + zetafun((2*k)-1,s)); %KHACK - SORTED
 PAB = 1-PAC; %KHACK - SORTED
 PCE = zetafun((2*k)-1,s)/( zetafun((2*k),s) + zetafun((2*k)-1,s)); %KHACK - SORTED
 PCF = 1-PCE; %KHACK - SORTED

%% Test form D

for n=1:7
    %%% these equations are valid for PBD but not PDB
    
   start(n)= floor(2^(n-20))+floor(2^(n-19))+floor(2^(n-18))...
        +floor(2^(n-17))+floor(2^(n-16))+floor(2^(n-15))...
        +floor(2^(n-14))+floor(2^(n-19))+floor(2^(n-13))...
        +floor(2^(n-12))+floor(2^(n-11))+floor(2^(n-10))...
        +floor(2^(n-9))+floor(2^(n-8))+floor(2^(n-7))...
        +floor(2^(n-6))+floor(2^(n-5))+floor(2^(n-4))...
        +floor(2^(n-3))+floor(2^(n-2))+floor(2^(n-1))...
        +floor(2^n)+(2*(k-1));
    
    stopnum(n)=floor(2^(n+1))+floor(2^(n-20))+floor(2^(n-19))+floor(2^(n-18))...
        +floor(2^(n-17))+floor(2^(n-16))+floor(2^(n-15))...
        +floor(2^(n-14))+floor(2^(n-19))+floor(2^(n-13))...
        +floor(2^(n-12))+floor(2^(n-11))+floor(2^(n-10))...
        +floor(2^(n-9))+floor(2^(n-8))+floor(2^(n-7))...
        +floor(2^(n-6))+floor(2^(n-5))+floor(2^(n-4))...
        +floor(2^(n-3))+floor(2^(n-2))+floor(2^(n-1))...
        +floor(2^n)-1+(2*(k-1));
   clear probstemp
    for m=start(n):stopnum(n)
        
        probstemp(m)=zetafun(m,s);
    end
    probsDG (n) = sum(probstemp);
    probsDB (n)= 1-(startnumerator+sum(probsDG(1:n))); %KHACK startnumerator
    transDG (n) = probsDG(n)/(probsDG(n)+probsDB(n));
end

 normprobsDG=probsDG/sum(probsDG);
 probTrans= transDG.*normprobsDG;
PDG=sum(probTrans);
PDB=1-PDG;

%% Probabilities from G

for n=k*2+1:2:stop %KHACK value for 3 needs to be changed for k parametes - SORTED
    
    probsGH(n-k*2)=zetafun(n,s);
    probsGI(n-k*2)=zetafun(n+1,s);
     
 
   transGH (n-k*2) = probsGH(n-k*2)/(probsGH(n-k*2)+probsGI(n-k*2));
   transGI(n-k*2)=1- transGH(n-k*2);
    
end

 normprobsGH=probsGH/sum(probsGH);
 probTransGH= transGH.*normprobsGH;
PGH=sum(probTransGH);
 
  normprobsGI=probsGI/sum(probsGI);
 probTransGI= transGI.*normprobsGI;
 PGI=sum(probTransGI);
 
 %% PBD
 
 clear prob
 clear normprob
 
 prob=zeros(stop,1);
normprob=zeros(stop,1);

for n=k*2+1:stop
    
    prob(n)=zetafun(n,s);
    normprob(n)=prob(n)/(1-(startdenom)); %KHACK this needs to be recalculated for K - SORTED 
    
    bivaltemp = (bival(n,(~isnan(bival(n,:)))));  %num2str(dec2bin(n));
    
  
    for i = (floor(size(bivaltemp, 2)/2))+k:(size(bivaltemp, 2)-k)
        onetranstemp(n) = onetranstemp(n) + (bivaltemp(i));
    end
    
    onetrans=onetrans+onetranstemp;
    
   
    
    zerotranstemp(n)=floor(size(bivaltemp, 2)/2)-k-onetranstemp(n);
    

    
end

probTransBDones= normprob.*onetranstemp;
Transfinalones=sum(probTransBDones);

probTransBDzeros= normprob.*zerotranstemp;
Transfinalzeros=sum(probTransBDzeros);

PBDones=Transfinalones/(Transfinalones+Transfinalzeros);

PBDzeros=1-PBDones;
    
PAB
PAC
PCE
PCF
PBDones
PBDzeros
PDB
PDG
PGH
PGI