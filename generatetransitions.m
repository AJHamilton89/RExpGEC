function finaloutputseq=generatetransitions(stages,k)
outputseq=[ ];
clear finaloutputseq

for n=0:stages-1
    sequence = generateExpGtrellis(2^n);
    outputseq=[outputseq ;sequence]; 
end

%define end states

finalstates=finalstateExpGtrellis(2^(stages-1));

% trellisplotfinal=digraph(finalstates(:,1)',finalstates(:,2)');
% h = plot(trellisplotfinal,'Layout','layered');
finalstates(:,1:2)=finalstates(:,1:2)+length(outputseq)-2^stages;

for m=1:length(outputseq)
outputseq(m,3)=m;
end

outputseq(:,4)=outputseq(:,3)+outputseq(:,1);

finaloutputseq=[outputseq(:,4),outputseq(:,3),outputseq(:,2)];

finaloutputseq=[finaloutputseq; finalstates];

leaves=setdiff(finaloutputseq(:,2),finaloutputseq(:,1));

if k>0
for n=1:length(leaves)
    extension=buildgenericUECtree(k);
    extension(:,1:2)=extension(:,1:2)+max(max(finaloutputseq))-1;
    extension(1,1)=leaves(n);
    extension(2,1)=leaves(n);
    finaloutputseq=[finaloutputseq; extension];
end
end

leaves=setdiff(finaloutputseq(:,2),finaloutputseq(:,1));

for n=1:length(finaloutputseq)
    finaloutputseq(n,4) = ismember (finaloutputseq(n,2),leaves);
   finaloutputseq(n,5)=0;
end

%% idenitfy the final transition states
leaftrans=finaloutputseq(finaloutputseq(:,4)==1,:);

%% sort the trellis diagram

finaloutputseq = sortrows(finaloutputseq,1);
% HACK - potentially sorted
% finaloutputseq (2,3)=0;
% finaloutputseq (1,3)=1;

finaloutputseq = sortrows(finaloutputseq,4);

% %% set output states to match symbols
% finalstateloc=find(finaloutputseq(:,4)==1);
% 
% for n=1:length(finalstateloc)
%     
%     finaloutputseq(finalstateloc(n),2)=floor(length(finaloutputseq)/100)+100+n;
%     
% end

trellisplot=digraph(finaloutputseq(:,1)'+1,finaloutputseq(:,2)'+1);
figure (1)
h = plot(trellisplot,'Layout','layered');
end

