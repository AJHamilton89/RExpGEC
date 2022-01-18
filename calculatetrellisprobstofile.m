%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Generate RExpGEC trellis probabilities function

%   Inputs
%
% Trellis - defined trellis based on depth, k, coderate
% K
% Depth
% s
%
%   Outputs
%
% codeword - Output Codeword
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[probs]=calculatetrellisprobstofile(k,p1,depth,codingrate,maxcodes,mintrans,processes)


  % Deal with parallel processing on slurm
    process = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    if isnan(process)
        process = 0;
        processes = 1;
    else
        if process >= processes
            error('process >= processes');
        end
    end


minimumtransitions = 0;

s=zeta_p1_to_s(p1);

num_test_symbols=100;

trellis=generatetransitionstrellis(k,depth,codingrate);







% codeword=zeros(1,((width(trellis)-4)*length(input))); %initialise RExpGEC codeword as R
startstate=0;%ceil(rand*2)-1; %random start state - maybe useful to keep random for testing
fromstate=startstate;

NewCol = zeros(size (trellis,1),1); % intialise new column of trellis to start counting
trellis = [trellis NewCol]; %now append it.

while minimumtransitions<mintrans
    
    symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);
    
    
    
    reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);
    
    input = reorderedcodeword_probs;
    
    
    for n=1:length(input)
        
        bit=input(n);
        
        possiblestatesid = trellis(:,1) == fromstate; %create subset of transitions
        possiblestates = trellis(possiblestatesid,:);
        
        actualtransitionid = possiblestates(:,3) == bit; %navigate actual state
        actualtransition = possiblestates(actualtransitionid,:);
        
        tostate=actualtransition(2);
        
        
        %     codeword(n*(width(trellis)-4)-(width(trellis)-5):n*(width(trellis)-4))=actualtransition(5:end); %populate RExpGEC - this should be paramaterisable
        %
        
        fromstate=tostate; %create new from state to go back through loop
        
        [tf, index]=ismember(actualtransition,trellis,'rows'); %find the transition in the trellis
        trellis(index,width(trellis))=trellis(index,width(trellis)) + 1; %+1 on the transition
        
    end
    
    minimumtransitions=min(trellis(:,7));
    
    assert((0 <= tostate) && (tostate <= 1))
    
end
%% this section creates conditional probabilities
NewCol = zeros(size (trellis,1),1); % intialise new column of trellis to start counting
trellis = [trellis NewCol]; %now append it.

for n=1:2:size(trellis,1)
    
    trellis (n,size(trellis,2)) = trellis (n,(size(trellis,2)-1)) / ( trellis (n,(size(trellis,2)-1)) + trellis (n+1,(size(trellis,2)-1)) );
    trellis (n+1,size(trellis,2)) = 1 - trellis (n,size(trellis,2));
    
end

probs=trellis(:,width(trellis));





fn = sprintf('VariablesStorage/ProbsDepth=%i_Rate=%i_K=%i_maxcodes=%i_min_trans=%i_p1=%i.mat',depth,codingrate,k,maxcodes,mintrans,p1);
                save(fn,'probs')


end