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
clear all

% set REXpGEC parameters

for codingrate=2:4
    
    for depth=1:3
        
        for k=0:1:2
                                    
            for s = 1.1:0.1:3
                
                %variables that don't loop
                maxcodes = 1000; % set what the maximum value of codeset is. - Powers of 2  make sense
                num_symbols=200; %defines block size
                num_test_symbols=10000;
                
                %set EXIT parameters
                IA_count = 50; % Choose how many points to plot in the EXIT functions
                block_count = 10;
                
                
                % calculate dependent parameters
                
                
                %Generate the trellis
                trellis=generatetransitionstrellis(k,depth,codingrate);
                
                %% this section is for calculating the probabilities of the trellis - in practice this would be from an equation
                
                %Generate array of symbols in the zeta distribution
                symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);
                
                %Generate a reordered ExpG codeword from the symbols
                reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);
                
                
                %work through trellis and return how many transitions were recorded
                probs=calculatetrellisprobs(trellis,reorderedcodeword_probs);
                
                %initalise variables
                IAs = 0.999*(0:1/(IA_count-1):1);
                IE_av = zeros(1,IA_count);
                IE_means_hist = zeros(1,IA_count);
                area=0.0;
                
                % Determine each point in the EXIT functions
                for IA_index = 1:IA_count
                    
                    
                    
                    
                    
                    
                    
                    for block_index = 1:block_count
                        %Generate array of symbols in the zeta distribution
                        symbols=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);
                        
                        %Generate a reordered ExpG codeword from the symbols
                        clear reorderedcodeword
                        reorderedcodeword = generate_RExpGcodeword(k,symbols);
                        
                        %Generate RExpGEC codeword
                        clear RExpGEC
                        RExpGEC=generateRExpGEC(trellis,reorderedcodeword);
                        
                        %create apriori LLRs
                        clear apriori
                        apriori = generate_llrs(RExpGEC,IAs(IA_index));
                        
                        %create extrinsic LLRs
                        clear extrinsic
                        extrinsic = RExpGEC_trellis_decoder(apriori,trellis,probs); % this is measuring the MI of the output, i.e. the RExpG codewords, not the RExpGEC codeword
                        
                        %measure the Information
                        IEs_av(block_index) = measure_mutual_information_averaging(extrinsic);
                        IEs_hist(block_index) = measure_mutual_information_histogram(extrinsic,RExpGEC);
                    end
                    
                    
                    
                    
                    
                    %Store the mean and standard deviation of the results
                    IE_av(IA_index) = mean(IEs_av);
                    IE_means_hist(IA_index) = mean(IEs_hist);
                    
                    if(IA_index > 1)
                        area = area + (IE_av(IA_index)+IE_av(IA_index-1))*(IAs(IA_index)-IAs(IA_index-1))/2;
                    end
                    
                end
                
                %Create a figure to plot the results.
                figure;
                axis square;
                EXITtitle=sprintf('EXIT Function of RExpGEC  Depth=%i Rate=%i K=%i maxcodes=%i num symbols=%i s=%i.fig',depth,codingrate,k,maxcodes,num_symbols,s);
                title(EXITtitle);
                ylabel('I_A');
                xlabel('I_E');
                xlim([0,1]);
                ylim([0,1]);
                
                hold on;
                
                cnmean=IE_av;
                
                % % Plot the  EXIT function for the RExpgEC
                %         plot(IE_av,IAs,'-');
                plot(IAs,IE_means_hist,'r');
                plot(IAs,IE_av,'b');
                legend({'True quality','Claimed quality'},'Location','northwest');
                
                % Display the area beneath the EXIT function
                annotation('textbox','String',{['Area = ', num2str(area)]},'LineStyle','none','Position',[0.7 0.1 0.2 0.1]);
                
                
                hold on;
                
                fn1 = sprintf('Figures/EXITchartRExpGECDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_s=%i.fig',depth,codingrate,k,maxcodes,num_symbols,s);
                saveas(gcf,fn1)
                
                fn2 = sprintf('VariablesStorage/EXITRExpGECWorkspaceDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_s=%i.fig',depth,codingrate,k,maxcodes,num_symbols,s);
                save(fn2)
            end
            
        end
        
    end
    
end
