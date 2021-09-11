%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Generate RExpGEC codeword function   

%   Inputs
%
% Trellis - defined trellis based on depth, k, coderate
% Input - Array of bit values corresponding to input RExpG codeword
%
%   Outputs
%
% codeword - Output Codeword
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function codeword=generateRExpGEC(trellis,input)
clear codeword

codeword=zeros(1,((width(trellis)-4)*length(input))); %initialise RExpGEC codeword as R 
startstate=0;%ceil(rand*2)-1; %random start state - maybe useful to keep random for testing
fromstate=startstate;

for n=1:length(input)
    
    bit=input(n);
    
    possiblestatesid = trellis(:,1) == fromstate; %create subset of transitions
    possiblestates = trellis(possiblestatesid,:);
    
    actualtransitionid = possiblestates(:,3) == bit; %navigate actual state
    actualtransition = possiblestates(actualtransitionid,:);
     
    tostate=actualtransition(2);
    
    
    codeword(n*(width(trellis)-4)-(width(trellis)-5):n*(width(trellis)-4))=actualtransition(5:end); %populate RExpGEC - this should be paramaterisable
    
   
    fromstate=tostate; %create new from state to go back through loop
   

    
    
    
end

assert((0 <= tostate) && (tostate <= 1))