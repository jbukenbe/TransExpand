bus_loc = zeros(309,2);


for idx = 1:309
    this_bus = test(test(:,1)==idx,:);
    
    if ~isempty(this_bus)
        this_bus = this_bus(this_bus(:,2)>0,:);
        if ~isempty(this_bus)
            bus_loc(idx,1) = mean(this_bus(:,2));
            bus_loc(idx,2) = mean(this_bus(:,3));
        end
    end  
end