function plot2DBERfrom3Dcomplexity(x_mat,y_mat,z_mat,targetcomplexity)

complexity = y_mat(:,1);

[~,~,idx]=unique(round(abs(complexity-targetcomplexity)),'stable');
minVal=complexity(idx==1);



plot (y_mat(index,:),z_mat(index,:))


end

function plot2DBERfrom3DSNR(x_mat,y_mat,z_mat,targetSNR)


end