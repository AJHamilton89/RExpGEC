


function x = generate_random_symbols(a, distribution, parameter)

    x = zeros(1,a);
    random_number = rand(size(x));
    sum_prob = 0;
    symbol_value = 0;

    while sum(random_number > sum_prob) > 0
        symbol_value = symbol_value + 1;
        x(random_number > sum_prob) = symbol_value;
        sum_prob = sum_prob + feval(distribution, symbol_value, parameter);
    end
end

function Px = distribution_zeta(x, s)

    Px = x.^-s/zeta(s);
end



