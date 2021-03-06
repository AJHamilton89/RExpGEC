% EXIT function for a URC code used as an inner code
% Copyright (C) 2021  Alex Hamilton

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.



% Number of bits to encode
bit_count=100000;

% Number of a priori mutual informations to consider
IA_count=50;


% Channel SNR in dB
SNR = 2:0.2:2.2;


for SNRindex=1:length(SNR);
    
    EbN0=SNR(SNRindex)-10*log10(0.7669);
    % Noise variance
    N0 = 1/10^(SNR(SNRindex)/10);
    
    % Generate some random bits
    uncoded_bits  = round(rand(1,bit_count));
    coded_bits = URC4_encoder(uncoded_bits);
    
    
    % QPSK modulator
    % tx1 = -2*(coded_bits-0.5);
    
    %shape the bits
    Txbits_QPSK = reshape(coded_bits,[2,length(coded_bits)/2]);
    
    %QPSK Modulate
    tx1 = modulate(Txbits_QPSK);
    
    
    % Channel
    % -------
    
    % Uncorrelated Rayleigh fading channel
    channel = sqrt(1/2)*(randn(size(tx1)))+1i*randn(size(tx1));
    
    % AWGN channel
    %channel = ones(size(tx1));
    
    noise = sqrt(N0/2)*((randn(size(tx1)))+1i*randn(size(tx1))); %% this should be moved to the IA (inner-most) loop
    
    % Generate the received signal
    rx1 = tx1.*channel+noise;
    
    
    
    % QPSK demodulator
    % apriori_channel_llrs = (abs(rx1+1).^2-abs(rx1-1).^2)/N0;
    LLRs_QPSK_tilde = soft_demodulate(rx1, channel, N0);
    apriori_channel_llrs = reshape(LLRs_QPSK_tilde,size(coded_bits));
    
    % Plot the LLR histograms
    % display_llr_histograms([apriori_channel_llrs],[coded_bits]);
    
    % A priori mutual informations to consider
    IA = 0.999*(0:1/(IA_count-1):1);
    
    % Initialise results
    IE_av=zeros(1,IA_count);
    IE_hist=zeros(1,IA_count);
    area=0.0;
    
    % Consider each a priori mutual information
    for IA_index = 1:IA_count
        
        % Generate the a priori LLRs having the a priori mutual information considered
        apriori_uncoded_llrs = generate_llrs((abs(uncoded_bits-1)), IA(IA_index));
        
        % Do the BCJR
        extrinisic_llrs = URC4_decoder_bcjr(apriori_uncoded_llrs,apriori_channel_llrs);
        
        
        
        % Measure the mutual information of the extrinsic LLRs
        IE_hist(IA_index) = measure_mutual_information_histogram( extrinisic_llrs, uncoded_bits);
        IE_av(IA_index) = measure_mutual_information_averaging( extrinisic_llrs);
        
        % Update the area beneath the EXIT function
        if(IA_index > 1)
            area = area + (IE_av(IA_index)+IE_av(IA_index-1))*(IA(IA_index)-IA(IA_index-1))/2;
        end
    end
    
    
    % Plot EXIT function
    figure
    axis square
    xlim([0 1]);
    ylim([0 1]);
    ylabel('I_A');
    xlabel('I_E');
    title(['EXIT function for SNR = ', num2str(SNR(SNRindex)), ' dB,EbN0 = ', num2str(EbN0), ' dB']);
    hold on
    plot(IA,IE_hist,'r');
    plot(IA,IE_av,'b');
    legend({'True quality','Claimed quality'},'Location','northwest');
    
    % Display the area beneath the EXIT function
    annotation('textbox','String',{['Area = ', num2str(area)]},'LineStyle','none','Position',[0.7 0.1 0.2 0.1]);
    
    fn1 = sprintf('Figures/EXITURC_SNR=%i.fig',SNR(SNRindex));
    saveas(gcf,fn1)
    
    fn2 = sprintf('VariablesStorage/EXITURC4Workspace_SNR=%i.mat',SNR(SNRindex));
    save(fn2)
    
end