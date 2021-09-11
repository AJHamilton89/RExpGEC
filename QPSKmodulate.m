% bits is a 1xk vector of bits
% tx is a complex symbol
function txoutput = QPSKmodulate(totalbits)

txoutput = [];

% Specify the constellation points and the bit mapping here
constellation_points = [+1+1i; -1+1i; -1-1i; +1-1i]/sqrt(2);
bit_labels = [0,0; 0,1; 1,1; 1,0];

% Determine the number of bits per symbol and the number of constellation points here
k = size(bit_labels,2);
M = 2^k;

for n=1:k:length(totalbits)
    
    bits = totalbits(n:n+(k-1))';
    
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
    
    txoutput = [txoutput tx];
    
end
end
