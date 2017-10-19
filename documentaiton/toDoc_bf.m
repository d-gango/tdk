% átalakítja a vágólapon található stringet és visszamásolja a vágólapra
function toDoc_bf()
    str1  = clipboard('paste');
    expression1 = 'bm';  %\text{} kiszedés
    rep1 = 'mathbf';
    string1 = regexprep(str1,expression1,rep1);

    clipboard('copy',string1)
end
