%% Simulated Annealing
%{
runs = 50;
time = 1000;

val_log = zeros(time,runs);
best_val_log = zeros(runs,1);


for j = 1:runs
solut = round(rand(24,1));
temp = ones(time,1);
cur_val = valcalc(solut);
val_log(1,j) = cur_val;
best_val = 0;
best_sol = solut;

for i = 1:(time-1)
    temp(i:time) = temp(i:time).*.99;
    
    if mod(i,100) > 70
        temp(i:time) = temp(i:time).*1.03;
    end
    
    idx = randi(24);
    candsol = solut;
    candsol(idx) = abs(candsol(idx)-1);
    
    cand_val = valcalc(candsol);
    
    if best_val < cand_val
        best_sol = candsol;
    end
    
        
    if cur_val < cand_val 
        solut = candsol;
        cur_val = cand_val;
        best_val = max(best_val, cand_val);
        val_log(i+1,j) = cand_val;
    elseif rand() < temp(i)
        solut = candsol;
        cur_val = cand_val;
        val_log(i+1,j) = cand_val;
    else
        val_log(i+1,j) = cur_val;
    end
        
end
best_val_log(j) = best_val;
plot(val_log(:,j))
hold on
end
yyaxis right
plot(temp)
mean(best_val_log)
%}








%% Genetic Algo

%{
gen_lim = 500;
pop_size = 500;
migrant_size = 10;
gen_avg = zeros(gen_lim,1);
gen_max = zeros(gen_lim,1);
gen_min = zeros(gen_lim,1);
gen_diversity = zeros(gen_lim,1);
best_sol = zeros(1,54);
best_val = 0;
count = 0;

pop = round(rand(pop_size,54));
pop_vals = zeros(pop_size,1);

for v_idx = 1:pop_size
    pop_vals(v_idx) = valcalc(pop(v_idx,:));
end

[best_val,max_id] = max(pop_vals);
best_sol = pop(max_id,:);

for g_idx = 1:gen_lim   
    pop_samp_pi =  (1+pop_vals-min(pop_vals))./max(pop_vals);
    pop_samp_pi(1+pop_size-migrant_size) = max(pop_samp_pi)*.9;
    sumpi = sum(pop_samp_pi);
    cumsumpi = cumsum(pop_samp_pi./sumpi);

    new_pop = zeros(pop_size,54);
    for p_idx = 1:(pop_size-migrant_size)
        mum = find(rand()<cumsumpi,1);
        dad = find(rand()<cumsumpi,1);
        genome = [pop(mum,:) ; pop(dad,:)];

        new_pop(p_idx,:) = diag(genome(ceil(rand(1,54)*2),1:54));
        
        if rand()<.05
            mutate_id = randi(54);
            new_pop(p_idx,mutate_id) = 1-new_pop(p_idx,mutate_id);
        end
    end
    new_pop((1+pop_size-migrant_size):pop_size,:) = round(rand(migrant_size,54));
       
    pop = new_pop;
    
    
    for v_idx = 1:pop_size
        pop_vals(v_idx) = valcalc(pop(v_idx,:));
    end
    count = count + sum(pop_vals == 617);
    
    if max(pop_vals) > best_val
        [best_val,max_id] = max(pop_vals);
        best_sol = pop(max_id,:);
    end   
    
    gen_avg(g_idx) = mean(pop_vals);
    gen_max(g_idx) = max(pop_vals);
    gen_min(g_idx) = min(pop_vals);
    covnorm = ((pop)*(pop)')./24;
    gen_diversity(g_idx) = mean(mean(covnorm./diag(covnorm)')); 
end
plot(gen_max);
hold on
plot(gen_min);
plot(gen_avg);
yyaxis right
plot(gen_diversity)
best_val
count
%}




%% Not dumb search

sol = round(rand(1,54));
local_idx = 0;
best_val = valcalc(sol);
val_log = zeros(100,1);
calcs = 0;

for g_idx = 1:100
    
    
    if local_idx == 0
        step_vals = zeros(54,1);
        new_step = abs(eye(54) - sol);  
        for v_idx = 1:54
            step_vals(v_idx) = valcalc(new_step(v_idx,:));
            calcs = calcs +1;
        end
        
    elseif local_idx == 1
        step_vals = zeros(54^2,1);
        new_step = abs(min(ones(54^2,54),repelem(eye(54),54,1)+repmat(eye(54),54,1))-sol);
        for v_idx = 1:length(new_step) 
            step_vals(v_idx) = valcalc(new_step(v_idx,:));
            calcs = calcs +1;
        end 
    end
    
    
    
    
    if max(step_vals) > best_val    
        [best_val, best_id] = max(step_vals);
        sol = new_step(best_id,:);
        local_idx = 0;
    elseif best_val < 617
        local_idx = local_idx+1
    end

    val_log(g_idx) = best_val;
end
plot(val_log)
calcs






function v = valcalc(x)
costm = [10	10	10	-5	-5	-5	-40
9	8	7	-6	-7	-8	20
5	5	5	5	5	5	5
1	2	3	-20	-20	-20	100
-6	-6	-6	8	8	8	-8
-8	-3	5	10	4	8	-30
3	4	4	0	0	-8	20
20	1	1	-30	-30	-30	60
8	4	1	-3	-4	5	10
20	-8	-15	5	8	-5	10
13	12	-20	-5	-8	-8	14
1	8	1	-7	4	5	-16
-5	-5	10	-5	0	-20	20
20	-20	5	-6	-30	8	-40
80	30	60	5	-4	-20	20
70	50	32	-20	0	8	-40
-50	-60	40	8	-30	-30	80
10	-70	20	10	-3	5	40];




%{
costm = [10	10	10	-5	-5	-5	-40
9	8	7	-6	-7	-8	20
5	5	5	5	5	5	5
1	2	3	-20	-20	-20	100
-6	-6	-6	8	8	8	-8
-8	-3	5	10	4	8	-30
3	4	4	0	0	-8	20
20	1	1	-30	-30	-30	60];
%}

des = reshape(x,18,3);
des = [des, des(:,1).*des(:,2), des(:,1).*des(:,3), des(:,2).*des(:,3), des(:,1).*des(:,2).*des(:,3)];

v = sum(sum(des.*costm));


end