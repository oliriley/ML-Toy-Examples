close all
clc

%need to load a data file from the stim optimization task
%expects a matrix of the following column, arbitrary rows
%           [1]             [2]    [3] [4] [5] [6]        [7]       [8]
% Columns: Time(msec)---Successes---x---y---z---x-force---y-force---z-force

pathString = 'C:\Users\Oli\Google Drive\SPR 16\Capstone Data';
loadString = '122_20160427.txt';

addpath(pathString);
data = importdata(loadString,',');
data = data.data(2:end,:);
rmpath(pathString);

%time = data(:,1)/1000;
%figure(1)
%plot(time,data(:,2))

% figure(2)
% subplot(311)
% plot(time,data(:,3))
% subplot(312)
% plot(time,data(:,4))
% subplot(313)
% plot(time,data(:,5))

%color_step = ceil(length(data(:,1))/64);
%color_blks = 1:color_step:length(data(:,1))-color_step;
%map1 = colormap;

%figure(3)
%for colors = color_blks
    %plot3(data(colors:colors+color_step,3),data(colors:colors+color_step,4),data(colors:colors+color_step,5),'color',map1(find(color_blks == colors),:),'linewidth',2)
%    xlim([min(data(:,3)),max(data(:,3))])
%    ylim([min(data(:,4)),max(data(:,4))])
%    zlim([min(data(:,5)),max(data(:,5))])

%    drawnow
%    pause(0.01)
%    hold on
%end
%xlabel('Z')
%ylabel('X')
%zlabel('Y')
%bar = colorbar;
%caxis([0 max(time)])
%ylabel(bar,'Time (s)')
%hold off

dists = zeros(1,length(data(:,1)));
for x = 1:length(dists)
    dists(x) = sqrt(data(x,3)^2+data(x,4)^2+data(x,5)^2);
end
max(dists)
%%%%%%%%%%%%%%%
% %Experiment: cycle through colors for each point
% %Number of colors to cycle through
% numpoints = 64;
% % Make a colormap 
% map = [ zeros(numpoints/8,1), zeros(numpoints/8,1), linspace(0.5+8/numpoints,1,numpoints/8)'
%         zeros(numpoints/8,1), linspace(0+8/numpoints,.5,numpoints/8)', ones(numpoints/8,1)
%         zeros(numpoints/8,1), linspace(0.5+8/numpoints,1,numpoints/8)', ones(numpoints/8,1)
%         linspace(0+8/numpoints,.5,numpoints/8)', ones(numpoints/8,1), linspace(1-8/numpoints,.5,numpoints/8)'
%         linspace(0.5+8/numpoints,1,numpoints/8)', ones(numpoints/8,1), linspace(.5-8/numpoints,0,numpoints/8)'
%         ones(numpoints/8,1), linspace(1-8/numpoints,.5,numpoints/8)', zeros(numpoints/8,1)
%         ones(numpoints/8,1), linspace(.5-8/numpoints,0,numpoints/8)', zeros(numpoints/8,1)
%         linspace(1-8/numpoints,.5,numpoints/8)', zeros(numpoints/8,1), zeros(numpoints/8,1)];
% figure(4)
% plot3(0,0,0,'w')
% xlim([min(data(:,3)),max(data(:,3))])
% ylim([min(data(:,4)),max(data(:,4))])
% zlim([min(data(:,5)),max(data(:,5))])
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% hold on
% 
% for index = 1:numpoints-1:length(data(:,1))-numpoints
%     plot3(data(index:index+numpoints-1,3),data(index:index+numpoints-1,4),data(index:index+numpoints-1,5),'color',map(1 + mod(index,length(map(:,1))),:),'Linewidth',2)
%     drawnow
% end
% colorbar;
% hold off