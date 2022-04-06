% function for creating 2D BER plots with a target SNR from a 3D
% plot

function [x,y] = plot2DBERfrom3DSNR(x_mat,y_mat,z_mat,targetSNR)

complexity = y_mat(:,1); 

SNR = x_mat(1,:);

 [~,index] = (min(abs(SNR-targetSNR))) ;
 
 x= complexity;
 y = z_mat(:,index);

% figure 
% 
% semilogy (complexity,z_mat(:,index))
% 
% 
% ylabel('BER')
% xlabel('Complexity')
% pltitle = ['Complexity vs BER at SNR of ', num2str(SNR(index)) ,'dB'];
% title(pltitle)

end