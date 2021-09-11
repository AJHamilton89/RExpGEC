


%% Inputs
% a - how many symbols to generate
% Px - finite probably distribution

function x = generate_random_symbols_Px(a, Px)

    x = zeros(1,a);
    random_number = rand(size(x));
    sum_prob = 0;
    symbol_value = 0;
    max=size(Px,2);

    while sum(random_number > sum_prob) > 0
        symbol_value = symbol_value + 1;
        if symbol_value > max %should this be less than
            break
        end
        x(random_number > sum_prob) = symbol_value;
        sum_prob = sum_prob + Px(symbol_value);
    end
end