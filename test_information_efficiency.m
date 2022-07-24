%%%code to calcuate information efficiency of codewords
clear all

index = 1;

kmax=3;
p1_count=20;


% p1 = (0:(p1_count-1))/(p1_count-1);
p1 = 0.00001+0.9999*(0:1/(p1_count-1):1);

for max=[1000];%10,100,1000,
    
    for p1_index = 1:p1_count
        
        %set S parameter - only required for zeta
        s=zeta_p1_to_s(p1(p1_index));
        
        s_storage(p1_index)=s;
        
        %create probabibility distribution
        for x=1:max
            Px(x) = zetafun(x,s);
        end
        
       
%         %geometric
%         for x=1:max
%         Px (x)= p1(p1_index)*(1-p1(p1_index)).^(x-1);
%         end
        
         %normalise Probability distribution
        Px=Px./sum(Px);
        
        
        %         codewords = createExpGcodes(max,k);
        %
        %         codelengths = sum(~isnan(codewords),2);
        
        for k=0:kmax
            
            %             %set how many symbols we need to go down (also includes the test length)
            %             num_symbols = 1000;
            
            codewords = createExpGcodes(max,k);
            
            codelengths = sum(~isnan(codewords),2);
            
            
            
            
            
            %% create RExp codewords
            
            %             %MONTE CARLO approach
            %             %generate symbols
            %             symbols=generate_zeta_symbols_finite_dict(num_symbols,max,s);
            %
            %
            %             codewordout=[];
            %
            %             for num=1:length(symbols)
            %                 codeword=idcodeword(symbols(num),codewords);
            %                 codewordout = [codewordout,codeword];
            %             end
            %             averagecodewordlength=length(codewordout)/length(symbols);
            %
            %
            
            
            %zeta
            
            
            
            
            averagecodewordlength=sum(Px.*codelengths');
            
            information=measure_information_probdist(Px);
            
            %simulation
            
            rate(k+1)=information/averagecodewordlength;
            
        end
        ratenew(index,1:k+1)=rate;
        index=index+1;
        
    end
    figure
    
    
    plot(p1,ratenew)
    
    ylabel('\eta_0')
    title(['Information Efficiency vs P1 of zeta distributions L =' num2str(max)])
    xlabel('P1')
    ylim([0 1])
    xlim([0 1])
    grid on
    legend("k = 0","k = 1","k = 2","k = 3","k = 4","k = 5")
    drawnow
    
    figure
    
    plot(s_storage,ratenew)
    
    ylabel('\eta_0')
    title(['Information Efficiency vs S of zeta distributions L =' num2str(max)])
    xlabel('S')
    ylim([0 1])
    xlim([1 13])
    grid on
    legend("k = 0","k = 1","k = 2","k = 3","k = 4","k = 5")

    drawnow
    
    
end