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


codingrate=2;
num_test_symbols=10000;
p1_count=10;
p1 = 0.00001+0.9999*(0:1/(p1_count-1):1);

%set EXIT parameters
IA_count = 10; % Choose how many points to plot in the EXIT functions
block_count = 1;


% set REXpGEC parameters
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 1000; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=200; %defines block size

for k=1:2;
    
    
    %Generate the trellis
    trellis=generatetransitionstrellis(k,depth,codingrate);
    for p1_index = 1:p1_count
        
        s=zeta_p1_to_s(p1(p1_index));
        s_storage(p1_index)=s;
        
        
        
        
        
        
        
        
        
        %% this section is for calculating the probabilities of the trellis - in practice this would be from an equation
        
        %Generate array of symbols in the zeta distribution
        symbols_probs=generate_zeta_symbols_finite_dict(num_test_symbols,maxcodes,s);
        
        %Generate a reordered ExpG codeword from the symbols
        reorderedcodeword_probs = generate_RExpGcodeword(k,symbols_probs);
        
        
        %work through trellis and return how many transitions were recorded
        probs=calculatetrellisprobs(trellis,reorderedcodeword_probs);
        
        %initalise variables
        IAs = (0:(IA_count-1))/(IA_count-1);
        IE_means = zeros(1,IA_count);
        IE_stds = zeros(1,IA_count);
        
        % Determine each point in the EXIT functions
        for IA_index = 1:IA_count
            
            IEs = zeros(block_count);
            
            
            
            
            
            for block_index = 1:block_count
                %Generate array of symbols in the zeta distribution
                symbols=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);
                
                %Generate a reordered ExpG codeword from the symbols
                reorderedcodeword = generate_RExpGcodeword(k,symbols);
                
                %Generate RExpGEC codeword
                RExpGEC=generateRExpGEC(trellis,reorderedcodeword);
                
                %create apriori LLRs
                apriori = generate_llrs(RExpGEC,IAs(IA_index));
                
                %create extrinsic LLRs
                extrinsic = RExpGEC_trellis_decoder(apriori,trellis,probs); % this is measuring the MI of the output, i.e. the RExpG codewords, not the RExpGEC codeword
                
                %measure the Information
                IEs(block_index) = measure_mutual_information_averaging(extrinsic);
                %IEs(block_index) = measure_mutual_information_histogram(extrinsic,RExpGEC);
            end
            
            
            
            
            
            %Store the mean and standard deviation of the results
            IE_means(IA_index) = mean(IEs);
            IE_stds(IA_index) = std(IEs);
        end
        
        %Create a figure to plot the results.
        figure;
        axis square;
        EXITtitle=sprintf('EXIT Function of CND with Bands Depth=%i Rate=%i K=%i maxcodes=%i num symbols=%i s=%i.fig',depth,codingrate,k,maxcodes,num_symbols,s);
        title(EXITtitle);
        ylabel('I_A');
        xlabel('I_E');
        xlim([0,1]);
        ylim([0,1]);
        
        hold on;
        
        cnmean=IE_means;
        
        % % Plot the  EXIT function for the RExpgEC
        plot(IE_means,IAs,'-');
        plot(IE_means+IE_stds,IAs,'--');
        plot(IE_means-IE_stds,IAs,'--');
        hold on;
        
        fn1 = sprintf('Figures/EXITchartDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_s=%i.fig',depth,codingrate,k,maxcodes,num_symbols,s);
        saveas(gcf,fn1)
        
        
    end
end
