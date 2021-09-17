% QPSK modulator using natural mapping
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

% bits is a 1xk vector of bits
% tx is a complex symbol
function tx = modulate(bits)

    % Specify the constellation points and the bit mapping here
    constellation_points = [+1+1i; -1+1i; -1-1i; +1-1i]/sqrt(2);
    bit_labels = [0,0; 0,1; 1,1; 1,0];
    
    % Determine the number of bits per symbol and the number of constellation points here
    k = size(bit_labels,2);
    M = 2^k;
    
    % Check that all the vectors and matrices have the correct dimensions
    if ~isequal(size(constellation_points),[M,1]) || ~isequal(size(bit_labels),[M,k]) || size(bits,1)~=k
        error('wrong dimensions');
    end

    symbol_labels = zeros(M,1);
    symbols = zeros(1,size(bits,2));
    for bit_index = 1:k
        symbol_labels = symbol_labels + bit_labels(:,bit_index)*2^(bit_index-1);
        symbols = symbols + bits(bit_index,:)*2^(bit_index-1);
    end
    
    % Determine which symbol is mapped to the input bits
    tx = zeros(size(symbols));
    for symbol_index = 1:M
        tx(symbols == symbol_labels(symbol_index)) = constellation_points(symbol_index);
    end
end
    