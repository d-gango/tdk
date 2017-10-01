% átalakítja a vágólapon található stringet és visszamásolja a vágólapra
function toDoc()
    str1  = clipboard('paste');
    expression1 = '\x5Ctext\x7B(\w*)\x7D';  %\text{} kiszedés
    rep1 = '$1';
    string1 = regexprep(str1,expression1,rep1);
    
    expression2 = '(\w*)\x27';  % ' kiszedés
    rep2 = '\x5C\x64ot\x7B$1\x7D';
    string2 = regexprep(string1,expression2,rep2);

    old = {'(t)', 'l0', 'l1', 'l2', 'm0', 'm1', 'm2', 'I0', 'I1', 'I2', 'fi0', 'x0', 'y0', 'q0', 'q1', 'q2', '\dot{q0}', '\dot{q1}' ,'\dot{q2}'};
    new = {'', 'l_0', 'l_1', 'l_2', 'm_0', 'm_1', 'm_2', '\theta_0', '\theta_1', '\theta_2', '\varphi_0', 'x_0', 'y_0', 'q_0', 'q_1', 'q_2', '\dot{q}_0', '\dot{q}_1' ,'\dot{q}_2'};
    rep = replace(string2, old, new);
    clipboard('copy',rep)
end
