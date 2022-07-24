%check if EXIT Charts overlap
SNR = -5:0.1:5;

plotEXIT=1;
p1_count=10;
%p1 = 0.00001+0.9999*(0:1/(p1_count-1):1);
p1=0.1:0.1:0.9;

for num_symbols= [10000] %blocklength - just choose a large one
    
    for maxcodes = [100] %finite dictionary
        
        for depth = 1:2
            
            for codingrate = 2:3
                
                for k=1:2
                    
                    for p1index = 1:size(p1,2)
                        
                        
                        RExpGECfilename = sprintf('VariablesStorage/EXITRExpGECWorkspaceDepth=%i_Rate=%i_K=%i_maxcodes=%i_num_symbols=%i_p1=%i.mat',depth,codingrate,k,maxcodes,num_symbols,p1(p1index));
                        
                        load(RExpGECfilename,'IAs','IE_means_hist','area')
                        RExpGECx = IAs;
                        RExpGECy = IE_means_hist;
                        RExpGECarea=area;
                        
                        for URC = [2,4,8]
                            
                            for SNRindex = 1:size(SNR,2) %Rob suggests 7D matrix and break when
                                
                                
                                %load URC
                                URCfilename = sprintf('VariablesStorage/EXITURC%iWorkspace_SNR=%i.mat',URC,SNR(SNRindex));
                                load(URCfilename,'IA','IE_hist','area')
                                URCx=IE_hist;
                                URCy=IA;
                                URCarea=area;
                                
                                
                                %load RExpGEC
                                
                                
                                vq1 = interp1(RExpGECx, RExpGECy,sort([RExpGECx,URCx]));
                                % get the X axis from 'vanilla' of both and cross compare so do two
                                
                                % interpolations, each based on the true
                                vq2 = interp1(URCx, URCy, sort([RExpGECx,URCx]));
                                resultvector = vq1-vq2;
                                negativeelements = find(resultvector<0);
                                
                                Area=URCarea-(1-RExpGECarea);
                                
                                if sum(negativeelements) == 0
                                    %if it's the first open EXIT tunnel
                                    %break out and store the SNR value in
                                    %the array below
                                    
                                    break
                                end
                            end
                            
                            if plotEXIT ==1
                                
                                % Plot EXIT function
                                figure
                                axis square
                                xlim([0 1]);
                                ylim([0 1]);
                                ylabel('I_A');
                                xlabel('I_E');
                                title(['EXIT function for SNR = ', num2str(SNR(SNRindex)), ' dB, URC = ', num2str(URC) ' RExpGEC k =' num2str(k) ' p1 = ' num2str(p1(p1index))]);
                                hold on
                                plot([RExpGECy],[RExpGECx],'r');
                                plot([URCy],[URCx],'b');
                                
                                %Display the area beneath the EXIT function
                                annotation('textbox','String',{['Tunnel Area = ', num2str(Area)]},'LineStyle','none','Position',[0.6 0.1 0.2 0.1]);
                                
                                %Display the area beneath the EXIT function
                                annotation('textbox','String',{['Area under RExpgEC = ', num2str((1-RExpGECarea))]},'LineStyle','none','Position',[0.6 0.2 0.2 0.1]);
                                
                                s=zeta_p1_to_s(p1(p1index));
                                
                                [av_codeword_length,entropy] = Find_av_codeword_length_zeta(maxcodes,k,s);
                                
                                Efficiency = (entropy/av_codeword_length)/codingrate;
                                
                                %Display the area beneath the EXIT function
                                annotation('textbox','String',{['Effective R = ', num2str((Efficiency))]},'LineStyle','none','Position',[0.6 0.3 0.2 0.1]);
                                
                                
                            end
                            
                            EXITTunnellowestSNR(URC,depth,codingrate,k,maxcodes,num_symbols,p1index) = SNR(SNRindex); %check if this is right
                            
                        end
                    end
                end
            end
        end
    end
end
