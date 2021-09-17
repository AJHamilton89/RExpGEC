% Soft QPSK demodulator using natural mapping
% Copyright (C) 2010  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

% apriori_llrs is a 1xk vector of a priori LLRs
% rx is a complex symbol
% channel is a complex channel coefficient
% N0 is the noise power spectral density
% extrinsic_llrs is a 1xk vector of extrinsic LLRs
function extrinsic_llrs = soft_demodulate(rx, channel, N0)

    % Specify the constellation points and the bit mapping here
    constellation_points = [+1+1i; -1+1i; -1-1i; +1-1i]/sqrt(2);
    bit_labels = [0,0; 0,1; 1,1; 1,0];

    % Determine the number of bits per symbol and the number of constellation points here
    k = size(bit_labels,2);
    M = 2^k;
    N = length(rx);

    % Check that all the vectors and matrices have the correct dimensions
    if ~isequal(size(constellation_points),[M,1]) || ~isequal(size(bit_labels),[M,k])
        error('wrong dimensions');
    end

    % Calculate the log probability for each symbol
    symbol_probabilities = -abs(repmat(rx,M,1)-constellation_points*channel).^2/N0;

    % Calculate the extrinsic LLRs
    log_p0 = zeros(k,N);
    log_p1 = zeros(k,N);
    for bit_index = 1:k
        log_p0(bit_index,:) = maxstar(symbol_probabilities(bit_labels(:,bit_index) == 0,:));
        log_p1(bit_index,:) = maxstar(symbol_probabilities(bit_labels(:,bit_index) == 1,:));
    end
    extrinsic_llrs = log_p1-log_p0;




end
