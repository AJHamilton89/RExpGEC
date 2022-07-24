% Copyright (C) 2013  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

 
 
 
 
 
% Encoder function for a terminated unity-rate recursive convolutional code
% having 3 memory elements, a generator polynomial of [1,1,1] and a feedback
% polynomial of [1,1,0].

% As specified in <VLC>.P354

% b is a vector of uncoded bits
% c is a matrix of encoded bits
function c = URC8_encoder(b)

    % Initialise our output bit matrix
    c = zeros(1,length(b)); % One row because this is a unity-rate CC
    
    % We start in the all-zeros state
    s1 = 0;
    s2 = 0;
    s3 = 0;
    
    % Encode the uncoded bit sequence
    for bit_index = 1:length(b)
        
        % Determine the next state
        s1_plus = mod(b(bit_index)+s1+s2+s3, 2); % This uses the feedback polynomial
        s2_plus = s1;
        s3_plus = s2;


        % Determine the encoded bits
        c(1, bit_index) = mod(s1_plus, 2); % This uses the first generator polynomial
        
        % Enter the next state
        s1 = s1_plus;
        s2 = s2_plus;
        s3 = s3_plus;
    end
    
end