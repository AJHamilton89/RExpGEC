% %27/2/19 page 105 of logbook

% Alterations - 7/5/19

% Finally created empircal vs predicted results : caveat - only for k = 0,
% to speed up run time suggest store several large codeword sets
% (10000000+ codes) at varying values of s;

 clear all
depth = 2;
maxcodes = 32; % set what the maximum value of codeset is. - Powers of 2  make sense
k=0;
codewordout=[];
num_symbols=100000;
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

Filename = sprintf('codeword_%s_S=%d_size=%d.mat', datestr(now,'mm-dd-yyyy HH-MM'), s, num_symbols);
save(Filename,'codewordout')
%% section to calculate predicted results - assumption is to work without truncation first

maximum_depth=ceil(log2(maxcodes)); %work out what the true maximum depth is

% final_branches= [(maximum_depth-1)^2-1:maximum_depth^2-1]; MATHS was
% wrong...
trellis = NaN(maximum_depth,2^(maximum_depth-1));

for n=1:maximum_depth
temp_final_branches= [2^(n-1):2^(n)-1];
trellis(n,1:length(temp_final_branches)) = temp_final_branches;


end

% mother=1;
% [UECl,UECr,FLCl,FLCr]=findtrellisprobs(mother,trellis,Px)
% UEC trans is mother vs children
% FLC trans is sibling rivalry


% % mother = 1;
% % children=findchildren(mother,trellis);
% % UECtransterminal= sum(Px(mother))/(sum(Px(mother))+sum(Px(children)));
% % UECtranscarry=1-UECtransterminal;
% % 
% % firstchildren=findchild(mother,trellis);
% % left_children=findchildren(firstchildren(1),trellis);
% % right_children=findchildren(firstchildren(2),trellis);
% % 
% % FLCtransleft = sum(Px(left_children))/(sum(Px(left_children))+sum(Px(right_children)));
% % FLCtransright = 1-FLCtransleft;
% % 
% % 
% %                     % calculate first transition
% %                     left_symbols = [];
% %                     right_symbols = [];
% %                     for n=2:maximum_depth
% %                         
% %                     temp_final_branches = trellis(n,:);
% %                     temp_final_branches= temp_final_branches(~isnan(temp_final_branches));
% %                     lx = (length(temp_final_branches));
% %                     half = ceil(lx/2);
% % 
% %                     lower_half = temp_final_branches(1:half);
% %                     upper_half = temp_final_branches(half+1:end);
% %                     left_symbols = [ left_symbols, lower_half ];
% %                     right_symbols = [ right_symbols, upper_half ]; 
% %                     end
% % 
% %                     left_prob = sum(Px(left_symbols));
% %                     right_prob = sum(Px(right_symbols));
% % 
% %                     Conditional_left = left_prob/(left_prob+right_prob); % Value for first left turn
% %                     Conditional_right = right_prob/(left_prob+right_prob); % Value for first right turn
% % 
% %                     
                     % calculate second transition
                    
%% calculate transitions


Initial_state = 0; % set start state of trellis
current_state = Initial_state;

 ind0 = treematrix(:,4) == 1; % find leaf nodes
 num_leafs=sum(ind0);% how many leafs
%  reset_matrix=zeros(num_leafs,size(treematrix,2)); % create a matrix to add reset states
%  reset_matrix(:,1) = treematrix(ind0,2); %populate it
%  treematrix=[reset_matrix;treematrix]; %concatenate with actual matrix
 codewordcount=0;
% SORTED - note - causing errors as no return state!!!
treematrix=sortrows(treematrix,1);

for n=1:length(codewordout) % go through codeword input
  ind1 = treematrix(:,1) == current_state; % set index of row that is equal to current state
  ind2 = treematrix(:,3) == codewordout(n); % find index of rows with transition being correct
  ind3 = logical(ind2.*ind1); % find where both ind1 and ind2 are true
  actual_row = treematrix(ind3,:); % define actual rows  
  
  treematrix (ind3,5) =  treematrix (ind3,5) +1;   
  next_state = actual_row(:,2);
  if actual_row(:,4) == 1; % set condition to automatically return to 0 if end of codeword
      next_state = 0;
      codewordcount=codewordcount+1;
  end
 current_state=next_state;
      
end

%calculate empirical probabilities against total set
for n=1:length(treematrix)
    treematrix(n,6)=treematrix(n,5)/sum(treematrix(:,5));
end

treematrix(n,7)=0;
%calculate empirical probabilities against each transition
for n= 0: max(treematrix(:,2))
    ind1 = treematrix(:,1)==n;% find out which rows share common start node
    if any(ind1)
    comparison_rows = treematrix(ind1,:);
    
    comparison_rows(1,7)=comparison_rows(1,6)/(comparison_rows(1,6)+comparison_rows(2,6));
    comparison_rows(2,7)=1-comparison_rows(1,7);
    treematrix(ind1,:)=comparison_rows;
    end
end




if k==0 %currently only works for k=0


%put predicted results into treematrix
treematrix(:,8)=NaN;

%set column 8 to represent the true final nodes
value=1;

for n=1:length(treematrix)
if treematrix(n,4)==1
treematrix(n,8)=value;
value=value+1;
end
end

for n=1:length(treematrix) % calculate transistion probabilities for all except final probabilities - unless truncated
if ~isnan(treematrix(n,8))
mother=treematrix(n,8);
[UECl,UECr,FLCl,FLCr]=findtrellisprobs(mother,trellis,Px);
treematrix(n,9)=UECl;
treematrix(n-1,9)=UECr;

FLCs=find(treematrix(:,1)==treematrix(n-1,2));

treematrix(FLCs(1),9)=FLCl;
treematrix(FLCs(2),9)=FLCr;


end
end

%calculate differences
treematrix(:,10)=NaN;
treematrix(:,10)=abs(treematrix(:,9)-treematrix(:,7)).^2; %Squared error

treematrix(find(isnan(treematrix)))=0;
end

treematrix;


% Input state = 0
% If transition value = 1, next state = (look up input State = X, transition value =1)
% 
% If transition value = 0, next state = (look up input State = X, transition value =0)
% 
% Calculate Transitions (input state, output state) = Calculate Transitions (input state, output state) +1
% 
% if output state = = only child // BREAK
% 
% % set depth (calcluate end states/finite symbols)
% symbolnum=2^(depth+1)*(k+1);


