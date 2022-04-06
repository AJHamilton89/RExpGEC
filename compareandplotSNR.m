function compareandplotSNR(maxk,num_symbols,p1,codingrate,maxcodes,maxdepth,SNR)

figure
hold on


for k = 0:maxk
    for d = 1:maxdepth

        [a,b,c]=consolidate(k,num_symbols,p1,codingrate,maxcodes,d);
      
        [x,y] = plot2DBERfrom3DSNR(a,b,c,SNR);
       titlestring = ['RExpGEC k = ',num2str(k), ' depth =' ,num2str(d)];
       semilogy(x,y,'DisplayName',titlestring )
        
    end
   
    
    [a,b,c]=ExpGCCURCconsolidate(k,num_symbols,p1,codingrate,maxcodes,1);
    [x,y] = plot2DBERfrom3DSNR(a,b,c,SNR);
    
   titlestring = ['ExpG URC CC k = ',num2str(k)];
       semilogy(x,y,'DisplayName',titlestring )
 
end
 


ylabel('BER')
xlabel('Complexity')
pltitle = ['Complexity vs BER at fixed SNR of ', num2str(SNR), 'dB' ];
title(pltitle)
set(gca, 'YScale', 'log')
legend
 ylim([0.0001 1])

 