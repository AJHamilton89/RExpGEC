function [probs] = zetafun(symbolnum,s)

probs = (symbolnum).^-s/zeta(s);
end