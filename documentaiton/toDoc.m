% átalakítja a vágólapon található stringet és visszamásolja a vágólapra
function toDoc()
    str  = clipboard('paste');
    expression = '\x5Ctext\x7B(\w*)\x7D';  %\text{} kiszedés
    repl = '$1';

    string = regexprep(str,expression,repl);

    old = {'(t)', 'l0', 'l1', 'l2', 'm0', 'm1', 'm2', 'I0', 'I1', 'I2', 'fi0', 'x0', 'y0', 'q0', 'q1', 'q2'};
    new = {'', 'l_0', 'l_1', 'l_2', 'm_0', 'm_1', 'm_2', 'I_0', 'I_1', 'I_2', '\varphi_0', 'x_0', 'y_0', 'q_0', 'q_1', 'q_2'};
    rep = replace(string, old, new);
    clipboard('copy',rep)
end
