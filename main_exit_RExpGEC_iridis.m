
% Script for drawing the EXIT chart of the RExpGEC
% Copyright (C) 2021  Alex Hamilton
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

function main_exit_RExpGEC_iridis(processes)

if ~exist('processes','var')
    processes = 1;
end

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




% set REXpGEC parameters

%variables that don't loop

prev='';

p1_count=processes;
p1 = 0.00001+0.9999*(0:1/(p1_count-1):1);
p1_index=process+1;

maxcodes = 100; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=10000; %defines block size
num_test_symbols=10000;

%set EXIT parameters
IA_count = 50; % Choose how many points to plot in the EXIT functions
block_count = 1;



for codingrate=2:4
    
    for depth=1:2
        
        for k=1:2
            
            %Generate the trellis
            trellis=generatetransitionstrellis(k,depth,codingrate);
            
            
            
                
                s=zeta_p1_to_s(p1(p1_index));
                s_storage(p1_index)=s;
                
                %Generate array of symbols in the zeta distribution
                symbolsarray=generate_zeta_symbols_finite_dict(num_symbols*block_count,maxcodes,s);
                
                
                
                
                
                %% this section is for calculating the probabilities of the trellis - in practice this would be from an equation
                
                %Generate array of symbols in the zeta distribution
                symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);
                
                %Generate a reordered ExpG codeword from the symbols
                reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);
                
                
                %work through trellis and return how many transitions were recorded
                probs=calculatetrellisprobs(trellis,reorderedcodeword_probs);
                
                %initalise variables
                IAs = (0.999*(0:1/(IA_count-1):1))+0.0001;
                IE_av = zeros(1,IA_count);
                IE_means_hist = zeros(1,IA_count);
                IE_means_hist_noprobs = zeros(1,IA_count);
                area=0.0;
                
                % Determine each point in the EXIT functions
                for IA_index = 1:IA_count
                    
                    
                    
                    
                    
                    
                    for block_index = 1:block_count
                        
                        symbols=symbolsarray(((block_index-1)*num_symbols+1):((block_index)*num_symbols));
                        
                        
                        %Generate a reordered ExpG codeword from the symbols
                        clear reorderedcodeword
                        reorderedcodeword = generate_RExpGcodeword(k,symbols);
                        
                        %Generate RExpGEC codeword
                        clear RExpGEC
                        RExpGEC=generateRExpGEC(trellis,reorderedcodeword);
                        
                        %create apriori LLRs
                        clear apriori
                        apriori =  generate_llrs(RExpGEC,IAs(IA_index)); %is this the wrong LLR definition?
                        
                        %create extrinsic LLRs
                        clear extrinsic
                        clear exrinsicRExpG
                        [extrinsic,extrinsicRExpG] = RExpGEC_trellis_decoder(apriori,trellis,probs); % this is measuring the MI of the output, i.e. the RExpG codewords, not the RExpGEC codeword
                        
                        [extrinsicnoprobs,extrinsicRExpGnoprobs] = RExpGEC_trellis_decoder(apriori,trellis); % this is measuring the MI of the output, i.e. the RExpG codewords, not the RExpGEC codeword
                        
                        IEs_hist_noprobs(block_index) = measure_mutual_information_histogram(extrinsicnoprobs,RExpGEC);
                        %measure the Information
                        IEs_av(block_index) = measure_mutual_information_averaging(extrinsic);
                        IEs_hist(block_index) = measure_mutual_information_histogram(extrinsic,RExpGEC);
                    end
                    
                    
                    
                    
                    
                    %Store the mean and standard deviation of the results
                    IE_av(IA_index) = mean(IEs_av);
                    IE_means_hist(IA_index) = mean(IEs_hist);
                    
                    IE_means_hist_noprobs(IA_index) = mean(IEs_hist_noprobs);
                    
                    if(IA_index > 1)
                        area = area + (IE_means_hist(IA_index)+IE_means_hist(IA_index-1))*(IAs(IA_index)-IAs(IA_index-1))/2;
                    end
                    
                end
                
                %                 %Create a figure to plot the results.
                %                 figure;
                %                 axis square;
                %                 EXITtitle=sprintf('EXIT Function of RExpGEC  Depth=%i Rate=%i K=%i maxcodes=%i num symbols=%i p1=%i.fig',depth,codingrate,k,maxcodes,num_symbols,p1(p1_index));
                %                 title(EXITtitle);
                %                 ylabel('I_A');
                %                 xlabel('I_E');
                %                 xlim([0,1]);
                %                 ylim([0,1]);
                %
                %                 hold on;
                %
                %                 cnmean=IE_av;
                %
                %                 % % Plot the  EXIT function for the RExpgEC
                %                 %         plot(IE_av,IAs,'-');
                %                 plot(IAs,IE_means_hist,'r');
                %                 plot(IAs,IE_av,'b');
                %                 legend({'True quality','Claimed quality'},'Location','northwest');
                %
                %                 % Display the area beneath the EXIT function
                %                 annotation('textbox','String',{['Area = ', num2str(area)]},'LineStyle','none','Position',[0.7 0.1 0.2 0.1]);
                %
                %
                %                 hold on;
                %
                %                 fn1 = sprintf('Figures/EXITchartRExpGECDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_p1=%i.fig',depth,codingrate,k,maxcodes,num_symbols,p1(p1_index));
                %                 saveas(gcf,fn1)
                
                
                %Create a figure to plot the results.
%                 figure;
%                 axis square;
%                 EXITtitle=sprintf('Inverted EXIT Function of RExpGEC  Depth=%i Rate=%i K=%i maxcodes=%i num symbols=%i p1=%i.fig',depth,codingrate,k,maxcodes,num_symbols,p1(p1_index));
%                 title(EXITtitle);
%                 ylabel('I_A');
%                 xlabel('I_E');
%                 xlim([0,1]);
%                 ylim([0,1]);
%                 
%                 hold on;
%                 
%                 cnmean=IE_av;
%                 
%                 % % Plot the  EXIT function for the RExpgEC
%                 %         plot(IE_av,IAs,'-');
%                 plot(IE_means_hist,IAs);
%                 plot(IE_means_hist_noprobs,IAs,'b--o');
%                 %                 plot(IE_av,IAs,'b');
%                 %                 legend({'True quality','Claimed quality'},'Location','northwest');
%                 
%                 current_label = [prev; sprintf('   Probs K=%i p1=%i',k,p1(p1_index)); sprintf('No Probs K=%i p1=%i',k,p1(p1_index)) ];
%                 legend(current_label,'Location','northwest');
%                 prev = current_label;
%                 drawnow
%                 Display the area beneath the EXIT function
                                annotation('textbox','String',{['Area = ', num2str(1-area)]},'LineStyle','none','Position',[0.7 0.1 0.2 0.1]);
                
                                fn3 = sprintf('Figures/InvertedEXITchartRExpGECDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_p1=%i.fig',depth,codingrate,k,maxcodes,num_symbols,p1(p1_index));
                                saveas(gcf,fn3)
                
                                fn2 = sprintf('VariablesStorage/EXITRExpGECWorkspaceDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_p1=%i.mat',depth,codingrate,k,maxcodes,num_symbols,p1(p1_index));
                                save(fn2)
            end
            
        end
        
    end
    
end


