%Plotting Probabilities based on Bayes' rule - I think I've done my math
%wrong

clear all
close all
clc

X = .1:.1:.5;

bar = X;

deltas = zeros(length(X),length(X),length(X));
colors = 'brkcm';
for super = 1:length(X)
    P_X = X(super);
    figure
    plot(bar,P_X*ones(1,length(X)),'k*-'),hold on
    for mid = 1:length(X)
        P_Y = P_X+X(mid);
        for sub = 1:length(X)
            P_bar = bar(sub);
            deltas(super,mid,sub) = (P_X - P_Y*P_bar)/(1-P_bar);
        end
    end
    for mid = 1:length(X)
        plot(bar,squeeze(deltas(super,mid,:)),colors(mid)),hold on
        title(strcat('X = ', num2str(P_X)))
        xlabel('P(bar)')
        ylabel('P(foo|~bar)')
        ylim([-0.2,0.5])
        xlim([0.1,0.5])
        legend('X',strcat('Y = ',num2str(P_X+X(1))),strcat('Y = ',num2str(P_X+X(2))),strcat('Y = ',num2str(P_X+X(3))),strcat('Y = ',num2str(P_X+X(4))),strcat('Y = ',num2str(P_X+X(5))),'Location','SouthWest')
    end
end
