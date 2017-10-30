%Practice for Gradient Descent on a least-squares cost function in 1D
clear all
close all
clc

%Create random 'training data' by selecting true underlying parameters
%(slope and y-intercept) and then generating points with deviations from
%that line

%Select true parameters
P1_true = rand - 0.5;
P2_true = rand - 0.5;

%Set the 'difficulty' of the dataset by setting how large the variations
%will be and how many points the algorithm gets to work with
Spread = .3;
Size = 100;

%Set the learning rate
Aleph = 2;

%Choose points in the space to create data for
X_vals = 1*rand(Size,1);

%Create the training set
Train_set = zeros(size(X_vals));
for index = 1:Size
    Train_set(index) = (rand-0.5)*Spread + P1_true*X_vals(index) + P2_true;
end

%Generate an initial hypothesis
P1_guess = rand - 0.5;
P2_guess = rand - 0.5;

%Iteratively improve the hypothesis
J_old = 10;
J_new = 1;

iter = 0;
while (abs(J_old - J_new) > J_old * 1e-10)
    %Calculate the old cost function
    J_old = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set).^2);

    %Plot the current state
    plot(X_vals,Train_set,'o'), hold on
    plot([min(X_vals),max(X_vals)],[P1_guess*min(X_vals) + P2_guess,P1_guess*max(X_vals) + P2_guess])
    legend('Data',sprintf('J = %d',J_old))
    drawnow
    %pause(.05)
    
    %Update the parameters
    ddP1 = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set).*X_vals);
    ddP2 = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set));
    P1_guess = P1_guess - Aleph*ddP1;
    P2_guess = P2_guess - Aleph*ddP2;
    
    J_new = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set).^2);
    
    iter = iter + 1;
    if J_new>J_old
        Aleph = Aleph/2;
    end
end

plot(X_vals,Train_set,'o')
plot([min(X_vals),max(X_vals)],[P1_guess*min(X_vals) + P2_guess,P1_guess*max(X_vals) + P2_guess],'r','Linewidth',3)
plot([min(X_vals),max(X_vals)],[P1_true*min(X_vals) + P2_true,P1_true*max(X_vals) + P2_true],'g','Linewidth',3)
legend('Data',sprintf('J = %d',J_new),'Actual')
pause(1)
%%
%The same thing, without derivatives this time!
clear all
close all

disp('Section 2')

%Create random 'training data' by selecting true underlying parameters
%(slope and y-intercept) and then generating points with deviations from
%that line

%Select true parameters
P1_true = rand - 0.5;
P2_true = rand - 0.5;

%Set the 'difficulty' of the dataset by setting how large the variations
%will be and how many points the algorithm gets to work with
Spread = .5;
Size = 200;

%Set the learning rate
Aleph = 10;
Bet = 0.3;

%Choose points in the space to create data for
X_vals = 1*rand(Size,1);

%Create the training set
Train_set = zeros(size(X_vals));
for index = 1:Size
    Train_set(index) = (rand-0.5)*Spread + P1_true*X_vals(index) + P2_true;
end

%Generate an initial hypothesis
P1_guess = rand - 0.5;
P2_guess = rand - 0.5;

%Iteratively improve the hypothesis
J_old = 10;
J_new = 1;

iter = 0;
while ((abs(J_old - J_new) > J_old * 1e-7) && iter < 100000)
    %Calculate the old cost function
    J_old = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set).^2);

    %Plot the current state
    plot(X_vals,Train_set,'o'), hold on
    plot([min(X_vals),max(X_vals)],[P1_guess*min(X_vals) + P2_guess,P1_guess*max(X_vals) + P2_guess])
    legend('Data',sprintf('J = %d',J_old))
    drawnow
    %pause(.05)
    
    %Update the parameters
    P1_up = P1_guess + Bet;
    P1_down = P1_guess - Bet;
    P2_up = P2_guess + Bet;
    P2_down = P2_guess - Bet;
    
    P1_temp = P1_guess;
    P2_temp = P2_guess;
    if sum((P1_up*X_vals + P2_guess - Train_set).^2) < sum((P1_down*X_vals + P2_guess - Train_set).^2)
        P1_temp = P1_up;
    elseif sum((P1_up*X_vals + P2_guess - Train_set).^2) > sum((P1_down*X_vals + P2_guess - Train_set).^2)
        P1_temp = P1_down;
    end
    
    if sum((P1_guess*X_vals + P2_up - Train_set).^2) < sum((P1_guess*X_vals + P2_down - Train_set).^2)
        P2_temp = P2_up;
    elseif sum((P1_down*X_vals + P2_up - Train_set).^2) > sum((P1_guess*X_vals + P2_down - Train_set).^2)
        P2_temp = P2_down;
    end
    
    P1_guess = P1_temp;
    P2_guess = P2_temp;
    
    J_new = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set).^2);
    
    if J_new>J_old
        Aleph = Aleph/2;
        Bet = Bet/2;
    end
    
    iter = iter + 1;
    if mod(iter,200)==0
        Aleph = Aleph*2;
        Bet = Bet*2;
    end
    if Bet>=1
        Bet = 0.99;
    end
end

hold off
plot(X_vals,Train_set,'o'), hold on
plot([min(X_vals),max(X_vals)],[P1_guess*min(X_vals) + P2_guess,P1_guess*max(X_vals) + P2_guess],'r','Linewidth',3)
plot([min(X_vals),max(X_vals)],[P1_true*min(X_vals) + P2_true,P1_true*max(X_vals) + P2_true],'g','Linewidth',3)
J_actual = 1/length(X_vals) * sum((P1_true*X_vals + P2_true - Train_set).^2);
legend('Data',sprintf('J = %d',J_new),sprintf('Actual J = %d',J_actual))

%%
%Fun! Now do it in 3D!
clear all
close all

%Create random 'training data' by selecting true underlying parameters
%(slope in X and Y, and Z-intercept) and then generating points with deviations from
%that line EDIT: THIS ISN'T HOW 3D LINES WORK, IGNORE

%Select true parameters
P1_true = rand - 0.5;
P2_true = rand - 0.5;
P3_true = rand - 0.5;

%Set the 'difficulty' of the dataset by setting how large the variations
%will be and how many points the algorithm gets to work with
Spread = .3;
Size = 100;

%Set the learning rate
Aleph = .07;

%Choose points in the space to create data for
X_vals = 1*rand(Size,1);
Y_vals = 1*rand(Size,1);

%Create the training set
Train_set = zeros(Size,1);
for index = 1:Size
    Train_set(index) = rand*Spread + P1_true*X_vals(index) + P2_true*Y_vals(index) + P3_true;
end

%Generate an initial hypothesis
P1_guess = rand - 0.5;
P1_orig = P1_guess;
P2_guess = rand - 0.5;
P2_orig = P2_guess;
P3_guess = rand - 0.5;
P3_orig = P3_guess;

%Iteratively improve the hypothesis
J_old = 10;
J_new = 1;

iter = 0;
while (abs(J_old - J_new) > J_old * 0.0001)
    %Calculate the old cost function
    J_old = 1/Size * sum((P1_guess*X_vals + P2_guess*Y_vals + P3_guess - Train_set).^2);

    %Plot the current state
    plot3(X_vals,Y_vals,Train_set,'o'), hold on
    plot3([min(X_vals),max(X_vals)],[min(Y_vals),max(Y_vals)],[P1_guess*min(X_vals) + P2_guess*min(Y_vals) + P3_guess, P1_guess*max(X_vals) + P2_guess*max(Y_vals) + P3_guess])
    legend('Data',sprintf('J = %d',J_old))
    drawnow
    %pause(.05)
    
    %Update the parameters (I think this is how it works in 3D?)
    ddP1 = 1/Size * sum((P1_guess*X_vals + P2_guess*Y_vals + P3_guess - Train_set).*X_vals);
    ddP2 = 1/Size * sum((P1_guess*X_vals + P2_guess*Y_vals + P3_guess - Train_set).*Y_vals);
    ddP3 = 1/Size * sum((P1_guess*X_vals + P2_guess*Y_vals + P3_guess - Train_set));
    %Trying to derive in 3D - stalls in the same way as above, which is maybe a good sign?
%     ddP1 = 1/(2*Size) * sum(P1_guess*X_vals*2    + (P2_guess*Y_vals).^2 + P3_guess^2 + Train_set.^2 + 2*(P2_guess*X_vals.*Y_vals          + P3_guess*X_vals          + P2_guess*P3_guess*Y_vals - X_vals.*Train_set          - P2_guess*Y_vals.*Train_set - P3_guess*Train_set));
%     ddP2 = 1/(2*Size) * sum((P1_guess*X_vals).^2 + P2_guess*Y_vals*2    + P3_guess^2 + Train_set.^2 + 2*(P1_guess*X_vals.*Y_vals          + P1_guess*P3_guess*X_vals + P3_guess*Y_vals          - P1_guess*X_vals.*Train_set - Y_vals.*Train_set          - P3_guess*Train_set));
%     ddP3 = 1/(2*Size) * sum((P1_guess*X_vals).^2 + (P2_guess*Y_vals).^2 + P3_guess   + Train_set.^2 + 2*(P1_guess*P2_guess*X_vals.*Y_vals + P1_guess*X_vals          + P2_guess*Y_vals          - P1_guess*X_vals.*Train_set - P2_guess*Y_vals.*Train_set - Train_set));
    P1_guess = P1_guess - Aleph*ddP1;
    P2_guess = P2_guess - Aleph*ddP2;
    P3_guess = P3_guess - Aleph*ddP3;
    
    J_new = 1/length(X_vals) * sum((P1_guess*X_vals + P2_guess - Train_set).^2);
    
    iter = iter + 1;
    %if J_new>J_old
    %    Aleph = Aleph/2;
    %end
end

plot3(X_vals,Y_vals,Train_set,'o')
plot3([min(X_vals),max(X_vals)],[min(Y_vals),max(Y_vals)],[P1_guess*min(X_vals) + P2_guess*min(Y_vals) + P3_guess, P1_guess*max(X_vals) + P2_guess*max(Y_vals) + P3_guess],'r','Linewidth',3)
legend('Data',sprintf('J = %d',J_new))