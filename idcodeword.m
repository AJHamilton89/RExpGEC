function code = idcodeword(x,codewords)

symbolnumber=x;


code=codewords(symbolnumber,1:end);

code(isnan(code)) = [];
end