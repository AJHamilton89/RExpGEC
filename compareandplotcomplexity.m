function compareandplotcomplexity(maxk,num_symbols,p1,codingrate,maxcodes,maxdepth,complexity)

figure
hold on


for k = 0:maxk
    for d = 1:maxdepth

        [a,b,c]=consolidate(k,num_symbols,p1,codingrate,maxcodes,d);
      
        [x,y] = plot2DBERfrom3Dcomplexity(a,b,c,complexity);
       
     titlestring = ['RExpGEC k = ',num2str(k), ' depth =' ,num2str(d)];
       semilogy(x,y,'DisplayName',titlestring )
        
    end
   
    
    [a,b,c]=ExpGCCURCconsolidate(k,num_symbols,p1,codingrate,maxcodes,1);
    [x,y] = plot2DBERfrom3Dcomplexity(a,b,c,complexity);
    
     titlestring = ['ExpG URC CC k = ',num2str(k)];
       semilogy(x,y,'DisplayName',titlestring )
 
end
 


ylabel('BER')
xlabel('SNR')
pltitle = ['SNR vs BER at fixed complexity of ', num2str(complexity), ' ACS operations, p1 = ', num2str(p1), 'blocklength = ', num2str(num_symbols) ];
title(pltitle)
set(gca, 'YScale', 'log')
legend ('Location','southwest')
grid on
 ylim([0.0001 1])

 