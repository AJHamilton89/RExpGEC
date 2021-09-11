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

% The geometric distribution
% x is the symbol value
% p1 is the parameter of the distribution
% Px is the probability corresponding to the specified symbol value
function Px = distribution_geometric(x, p1)

    % See Eq (1) in (Maunder et al., 2013)
    Px = p1*(1-p1).^(x-1);
end