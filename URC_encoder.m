% Encoder function for a terminated unity-rate recursive convolutional code
% having 3 memory elements, a generator polynomial of [1,1,0,1] and a feedback
% polynomial of [1,0,1,1]. This is as used in the UMTS turbo code, as specified 
% in ETSI TS 125 212 (search for it on Google if you like). For more 
% information, see Section 2.2.1 of Liang Li's nine-month report 
% (http://users.ecs.soton.ac.uk/rm/wp-content/liang_li_nine_month_report.pdf)
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

% b is a vector of uncoded bits
% c is a matrix of encoded bits
function c = URC_encoder(b)

    % Initialise our output bit vectors
    c = zeros(size(b));
    
    % We start in the all-zeros state
    s1 = 0;
    
    % Encode the uncoded bit sequence
    for bit_index = 1:length(b)
        
        % Determine the next state
        s1_plus = mod(b(bit_index)+s1, 2); % This uses the feedback polynomial
	
        % Determine the encoded bits
        c(bit_index) = mod(s1_plus, 2); % This uses the first generator polynomial
        
        % Enter the next state
        s1 = s1_plus;
    end
    
end