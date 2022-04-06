% Script that consolidates a bunch of results files into one.
% On IRIDIS, this can be run using the commands:
%         module load matlab
%         matlab -nodesktop -r "consolidate"
% The resultant file results.mat can be loaded into Matlab using the command:
%         load('-MAT', 'results.mat', 'consolidated_results')

% Copyright (C) 2022 Alex Hamilton

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.

% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

function [x_mat,y_mat,z_mat] = ExpGCCURCconsolidate(k,num_symbols,p1,codingrate,maxcodes,depth)


format shortG

first = true;

filestring = ['ExpGCCURCresults_k=',num2str(k),'_R=',num2str(codingrate),'_p1=',num2str(p1),'_num_sym=',num2str(num_symbols),'_maxcode=',num2str(maxcodes),'_depth=',num2str(depth)];

% Get a list of all files in the current directory.
filenames = dir;
% For each file in the current directory...
for file_index = 1:length(dir)
    % ...see if it is one of our results files.
    if ~isempty(strfind(filenames(file_index).name, filestring)) %if ~isempty(strfind(filenames(file_index).name, 'results_k=',num2str(k),'_R=',num2str(codingrate),'_p1=',num2str(p1),'_num_sym=',num2str(num_symbols)))
        % If so, load the results.
        load('-MAT', filenames(file_index).name, 'results' , 'complexity_array');
    
        % Consolidate the results into a single matrix.
        if first
            consolidated_results = results;
            first = false;
        else
            consolidated_results = [consolidated_results;results];
        end
    end
end

% Sort the consolidated results and display them.
consolidated_results = sortrows(consolidated_results);
results = consolidated_results;


SNR_Vec=consolidated_results(:,1);

[x_mat,y_mat]=meshgrid(SNR_Vec,complexity_array); 

z_mat=consolidated_results(:,6:end)';

% figure
% surf(x_mat,y_mat,z_mat)
% set(gca, 'zScale', 'log')
% set(gca, 'yScale', 'log')
% zlim([0.001 1])
% ylabel ('Complexity')
% xlabel ('SNR')
% zlabel ('BER')

% Save the consolidated results into a binary file. This avoids the loss of precision that is associated with ASCII files.
save(['ExpGCCURCconsolidated',filestring,'.mat'], 'results', 'complexity_array', '-MAT');

% % Plot the consolidated results. This will be ignored on Lyceum.
% semilogy(consolidated_results(:,1),consolidated_results(:,4));
% title('BPSK modulation in an AWGN channel');
% ylabel('BER');
% xlabel('SNR (in dB)');

% xaxis is repmat of SNR
% zaxis is repmat of complexity
% yaxis is repmat of BER

end