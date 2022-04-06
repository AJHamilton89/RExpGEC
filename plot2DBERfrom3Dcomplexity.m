% function for creating 2D BER plots with a target complexity from a 3D
% plot

function [x,y] = plot2DBERfrom3Dcomplexity(x_mat,y_mat,z_mat,targetcomplexity)

complexity = y_mat(:,1);

 [~,index] = (min(abs(complexity-targetcomplexity))) ;

x=x_mat(index,:);
 
y=z_mat(index,:);

% figure
% 
% semilogy (x_mat(index,:),z_mat(index,:))
% 
% ylabel('BER')
% xlabel('Complexity')
% pltitle = ['SNR vs BER at fixed complexity of ', num2str(complexity(index)), 'ACS operations' ];
% title(pltitle)

end

