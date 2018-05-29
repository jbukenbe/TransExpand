% this function takes the problem input from matlab and formats it to run
% the opf problem in gams and output the cost from each scenario








    guel = (@(s,v) strcat(s,strsplit(num2str(v))));
    subs.name='subs';
    subs.uels = guel('s',1);
    
    lines.name = 'lines';
    lines.uels = guel('line',y);
    
    filename = ['MtoG', '.gdx'];
    wgdx(filename,subs);
    system 'gams "OPF_file" lo=3'
    
    dummygen.name='scen_cost';
    dummygen.form='sparse';
    output = rgdx('results',dummygen);
    scen_op_cost = output.val;