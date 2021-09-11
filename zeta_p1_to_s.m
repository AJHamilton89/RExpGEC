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

% Obtains the parameter s from the probability of a zeta-distributed symbol having the value 1
% p1 is the probability of a zeta-distributed symbol having the value 1
% s is the parameter of the zeta distribution
function s = zeta_p1_to_s(p1)

s = zeros(1,length(p1));
for index = 1:length(p1)


    if p1(index) <= 0
        s = 1;
    elseif p1(index) >= 1
        s = inf;
    else

        top = 2;
        while 1/zeta(top) < p1(index)
            top = top + 1;
        end
        bottom = top-1;

        while top-bottom > 0.00001
            middle = (bottom+top)/2;
            if 1/zeta(middle) < p1(index)
                bottom = middle;
            else
                top = middle;
            end
        end
        s(index) = middle;
    end
end