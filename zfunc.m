function [probs] = zfunc(symbolnum,s)

probs = (symbolnum+1).^-s/zeta(s);
end