for f_idx = 1:25
    filename = sprintf('%s_%d','run',f_idx);
    load(filename);
    est(f_idx) = transmission_problem.tot_cost;
    
    
    for s_idx = 1:10
        transmission_problem = tproblem;
        transmission_problem.scen_cluster = (1:10);
        
        
        
        scen_lines = transmission_problem.lines(transmission_problem.lines(:,4) == transmission_problem.scen_cluster(s_idx),1);
        l_idx = (1:51)+(51*(s_idx-1));
        scen_line_bin = zeros(1,51);
        scen_line_bin(scen_lines) = 1;
        spy_lines(f_idx,l_idx) = scen_line_bin;
        scen_cost(s_idx) = sum(line_cost(scen_lines));
        if s_idx
            scen_cost(s_idx) = scen_cost(s_idx)*(.5/9);
        else
            scen_cost(s_idx) = .5*scen_cost(s_idx);
        end
    end
        %outdata(f_idx) = [.5*ones(24,1);.5*ones(24*(9),1)./(9)]'*sum(true_cost,2)./24 + sum(scen_cost);
        %outdata(f_idx) = [.5*ones(24,1);.5*ones(24*(9),1)./(9)]'*true_cost./24 + sum(scen_cost);
end