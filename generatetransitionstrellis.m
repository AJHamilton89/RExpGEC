function trellis=generatetransitionstrellis(k,stages,codingrate) %need to add hamming distance parameter from same state
outputseq=[ ];
clear finaloutputseq
% set the coding rate is specifically 1/R, rather than R



for n=0:stages-1
    sequence = generateExpGtrellis(2^n);
    outputseq=[outputseq ;sequence];
end

%define end states

finalstates=finalstateExpGtrellis(2^(stages-1));


finalstates(:,1:2)=finalstates(:,1:2)+length(outputseq)-2^stages;

for m=1:length(outputseq)
    outputseq(m,3)=m;
end

outputseq(:,4)=outputseq(:,3)+outputseq(:,1);

treematrix=[outputseq(:,4),outputseq(:,3),outputseq(:,2)];

treematrix=[treematrix; finalstates];

leaves=setdiff(treematrix(:,2),treematrix(:,1));

if k>0
    for n=1:length(leaves)
        extension=buildgenericUECtree(k);
        extension(:,1:2)=extension(:,1:2)+max(max(treematrix))-1;
        extension(1,1)=leaves(n);
        extension(2,1)=leaves(n);
        treematrix=[treematrix; extension];
    end
end

leaves=setdiff(treematrix(:,2),treematrix(:,1));

for n=1:length(treematrix)
    treematrix(n,4) = ismember (treematrix(n,2),leaves);
    
end




%% sort the trellis diagram

treematrix = sortrows(treematrix,1);

treematrix = sortrows(treematrix,4);


% trellisplot=digraph(treematrix(:,1)'+1,treematrix(:,2)'+1);
% figure (1)
% h = plot(trellisplot,'Layout','layered');

%% turn to recursive trellis state

%double the trellis - create mirrored trellis
trellis=zeros(2*length(treematrix),width(treematrix));

trellis(1:length(treematrix),1:width(treematrix))=treematrix; %populate

trellis(length(treematrix)+1:end,1:width(treematrix))=treematrix; %populate

trellis(1:length(treematrix),1:2)=trellis(1:length(treematrix),1:2)*2; %set the first part to be one side of trellis

trellis(length(treematrix)+1:end,1:2)=(trellis(length(treematrix)+1:end,1:2)*2)+1;% set the second part to be the other side of the trellis


%loop back to original state and switch between top and bottom trellis when
%bit value = 1
for n=1:length(treematrix)
    
    %loop back to original state
    if trellis(n,4) == 1
        
        trellis(n,2)=0;
        
    end
    
    %switch to other trellis if encoded bit is 1
    if trellis(n,3) == 1
        
        trellis(n,2) = trellis(n,2) + 1;
        
    end
    
    
end

%bit value=0
for n=length(treematrix)+1:length(trellis)
    
    %loop back to original state
    if trellis(n,4) == 1
        
        trellis(n,2)=1;
        
    end
    
    %switch to other trellis if encoded bit is 1
    if trellis(n,3) == 1
        
        trellis(n,2) =    trellis(n,2) - 1;
        
    end
    
end

%%generate random transistion for UEC codeword

for n=1:2:length(trellis)/2
    
    transition = ceil(rand*(2^codingrate))-1; % generate random number 0-coding rate
    transitionbin = bitget(transition, codingrate:-1:1); %turn to binary
    transitionbinopp = ~transitionbin; %get the oppositie for the opposite side
    
    
    trellis(n,5:(5+codingrate-1))=transitionbin; % start populating trellis transitions - original transition
    
     trellis(n+1,5:(5+codingrate-1))=transitionbinopp; %ensure that transitions from the same state are orthogonal
    
    trellis(n+(length(trellis)/2),5:(5+codingrate-1))=transitionbinopp; %ensure that transitions accross trellis are orthogonal
    
     trellis(n+1+(length(trellis)/2),5:(5+codingrate-1))=transitionbin; %ensure that transitions from the same state are orthogonal
  
end


%% plot some things
% trellisplot=digraph(trellis(:,1)'+1,trellis(:,2)'+1);
% figure (2)
% h = plot(trellisplot,'Layout','layered');



end

