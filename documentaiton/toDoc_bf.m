% �talak�tja a v�g�lapon tal�lhat� stringet �s visszam�solja a v�g�lapra
function toDoc_bf()
    str1  = clipboard('paste');
    expression1 = 'bm';  %\text{} kiszed�s
    rep1 = 'mathbf';
    string1 = regexprep(str1,expression1,rep1);

    clipboard('copy',string1)
end
