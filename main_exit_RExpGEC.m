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

% set parameters

k=1;
depth=1; %interesting point that the encoder complexity is affected also
maxcodes = 1000; % set what the maximum value of codeset is. - Powers of 2  make sense
num_symbols=200; %defines block size
s = 3;
codingrate=2;
num_its=100;
num_test_symbols=10000;
num_runs=10;
block_size = 1; %set the block size (not based on no.symbols)
IA_count = 10; % Choose how many points to plot in the EXIT functions
block_count = 1;
frame_count = 1;

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




IAs = (0:(IA_count-1))/(IA_count-1);
IE_means = zeros(1,IA_count);
IE_stds = zeros(1,IA_count);

% Determine each point in the EXIT functions
for IA_index = 1:IA_count
    
    IEs = zeros(1,frame_count);
    
    %     This runs the simulation long enough to produce smooth EXIT functions.
    for frame_index = 1:frame_count
        
        %         a = round(rand(block_size-1,block_count));
        %         b = [a;mod(sum(a,1),2)];
        
        reorderedcodeword=NaN(2*codingrate*num_symbols,block_size); %create an array of NaNs much bigger than maximum block size
        RExpGEC=NaN(4*codingrate*num_symbols,block_size); %create an array of NaNs much bigger than maximum block size
        apriori=NaN(4*codingrate*num_symbols,block_size);
        
        for block_index = 1:block_count
            %Generate array of symbols in the zeta distribution
            symbols(:,block_index)=generate_zeta_symbols_finite_dict(num_symbols,maxcodes,s);
            
            %Generate a reordered ExpG codeword from the symbols
            reorderedcodewordtemp = generate_RExpGcodeword(k,symbols(:,block_index));
%             reorderedcodeword(1:length(reorderedcodewordtemp),block_index) = reorderedcodewordtemp;
            
            %Generate RExpGEC codeword
%             idnotNaN=(~isnan(reorderedcodewordtemp));
            RExpGECtemp=generateRExpGEC(trellis,reorderedcodewordtemp);
            RExpGEC(1:length(RExpGECtemp),block_index)=RExpGECtemp;
            
            
            apriori = generate_llrs(RExpGECtemp,IAs(IA_index));
%             extrinsic = zeros(size(apriori));
            
            
            
            
            extrinsic = RExpGEC_trellis_decoder(apriori,trellis,probs); % this is measuring the MI of the output, i.e. the RExpG codewords, not the RExpGEC codeword
%         extrinsic(:,block_index) = RExpGEC_trellis_decoder(apriori,trellis,probs); % this is measuring the MI of the output, i.e. the RExpG codewords, not the RExpGEC codeword

        end
        
        %IEs(frame_index) = measure_mutual_information_averaging(extrinsic);
                
        
        IEs(frame_index) = measure_mutual_information_histogram(extrinsic,RExpGECtemp); %HACK
        
        
    end
    %Store the mean and standard deviation of the results
    IE_means(IA_index) = mean(IEs);
    IE_stds(IA_index) = std(IEs);
end

%Create a figure to plot the results.
figure;
axis square;
title('EXIT Function of CND with Bands');
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





