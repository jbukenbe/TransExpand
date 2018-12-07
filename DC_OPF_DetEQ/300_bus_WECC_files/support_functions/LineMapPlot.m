% Plot lines over a map of the western US




%% Get coordinates for buses
load('bus_loc.mat');
bus_n = size(busxy,1);

known_coords = busxy(busxy(:,1)>0,:);
yrange = [max(known_coords(:,1)), min(known_coords(:,1))];
xrange = [max(known_coords(:,2)), min(known_coords(:,2))];


busxy(busxy(:,1)==0,1) = mean(yrange);
busxy(busxy(:,2)==0,2) = mean(xrange);

%% Get x and y coordinates for lines
load('cand_line_bus.mat');
line_n = size(cand_line_bus,1);

%% Scale coordinates to plot size
busxy(:,1) = (-1950.*(busxy(:,1)- yrange(2))./(yrange(1)-yrange(2)))*.85+1750;
busxy(:,2) = (1950.*(busxy(:,2)- xrange(2))./(xrange(1)-xrange(2)))*.9+200;
coordx = 195.*coordx;
coordy = 195.*coordy;


%% Plot map
A = imread('WECC_shell.png');
image(A);
xlim([0 1950]);
ylim([0 1950]);
hold on


%% Plot buses
scatter(busxy(:,2),busxy(:,1), 'fill')
%scatter(known_coords(:,2),known_coords(:,1), 'fill')

%% Plot lines
for l_idx = [1 29 43 19 13 15 3 10 27 26]
    %plot(coordx(l_idx,:),coordy(l_idx,:),'k','LineWidth',5);
    plot([busxy(cand_line_bus(l_idx,1),2),busxy(cand_line_bus(l_idx,2),2)],...
        [busxy(cand_line_bus(l_idx,1),1),busxy(cand_line_bus(l_idx,2),1)]);
end

