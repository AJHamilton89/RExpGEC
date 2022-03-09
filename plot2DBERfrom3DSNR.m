% function for creating 2D BER plots with a target SNR from a 3D
% plot

function plot2DBERfrom3DSNR(x_mat,y_mat,z_mat,targetSNR)

complexity = y_mat(:,1); 

SNR = x_mat(1,:);

 [~,index] = (min(abs(SNR-targetSNR))) ;

figure 

semilogy (complexity,z_mat(:,index))


ylabel('BER')
xlabel('Complexity')
pltitle = ['Complexity vs BER at SNR of ', num2str(SNR(index)) ,'dB'];
title(pltitle)

end