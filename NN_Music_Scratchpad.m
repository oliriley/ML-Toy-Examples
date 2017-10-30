%@@TODO: Implement biases, and the training thereof
%@@TODO: Output weights over coarse digitization instead of single number?
%@@TODO: Does normalizing help? Have implemented it for now.

clear all
close all

% Set parameters
compressionFactor = 4; % Retain 1/cF of the song data
secondsBack = 10; % Length of time to be used as input
learnRate = 0.1; % Smaller = slower adjustment of weights
trainFrac = 0.8; % Portion of data to be devoted to training set
trainingCycles = 1000; % # of times training data is run through 
                    % Currently arbitrary - find by looking at error trace
errTrace = zeros(1,trainingCycles); %To store error on validataion set over time
spacingFactor = 2000;

%First grab the song
[song,Fs] = audioread('C:\Users\Oli\Downloads\Simple Gifts.mp3');
% Make it one-channel
song = mean(song,2);

% 'Simple Gitfts' has a chunk in front of zeros, discard them
song = song(576:end);

%Compress by picking out every (compressionFactor)-th data point
song = song(1:compressionFactor:end);
Fs = Fs/compressionFactor;

% Normalize the song on (-1,1)
maxSong = max(song);
minSong = min(song);

song = 2*(song+minSong)/(maxSong+minSong)-1;

% sound(song,Fs)
% clear sound

%Random initial weights on (-1,1)
%
% Given a number of seconds to be used for input,
% find the length of the input vector whose components will hop backward a
% number of datapoints equal to a power of two, in a number of hops equal
% to that power, starting at 2^-1.
% 
% That vector length defines the sizes of the weight matrices

powerof2 = 1;
inputLength = 1;
x = 1;
while (x(end) < Fs*secondsBack)
    for index = 1:powerof2
        x(length(x)+1) = x(end) + 2^(powerof2);
        inputLength = inputLength+1;
        if x(end) > Fs*secondsBack
            break
        end
    end
    powerof2 = powerof2 + 1;
end

layers = ceil(log2(inputLength));

weights{1} = 2*(rand(round(inputLength/2),inputLength))-1;
biases{1} = 2*(rand(length(x),1))-1;
tempWeights{1} = weights{1}*0;
for index = 2:layers
    weights{index} = 2*(rand(ceil(length(weights{index-1}(:,1))/2),ceil(length(weights{index-1}(1,:))/2)))-1;
    tempWeights{index} = weights{index}*0;
    
    biases{index} = 2*(rand(ceil(length(weights{index-1}(:,1))/2)))-1;
end

% w = weights{layers};
% for subindex = 1:layers-1
%     w = w*weights{layers-subindex};
% end

% Create an example song w/ the random weights
% Choose random chunk of song as a seed.
m = randi([x(end) + 1, length(song)]);
randSong = song(m-x(end) : m);
 
newSong = randSong; %This will be used to compare outputs later
PRETTY SURE THIS IS TOTES WRONG - CHECK SIZES OF BIASES{#}
for index = 1:x(end)*3 %how many song datapoints to generate
    prompt = nan(inputLength,1);
    for subindex = 1:inputLength
        prompt(inputLength-subindex+1) = randSong(end - x(subindex));
    end
    
    for subindex = 1:layers
        prompt = weights(layers-subindex+1)*prompt + biases(layers-subindex+1);
        prompt = 2./(1+exp(-1*prompt))-1;
    end
    
    randSong(length(randSong)+1) = prompt;
end

% Split data into training and validation datasets
trainData = [];
trainAns = [];

validData = [];
validAns = [];

% Assign datapoints randombly between training and validation, in
% proportion to trainFrac
for index = 1:round((length(song)-x(end))/spacingFactor)
    disp(strcat('Constructing datasets:`',num2str(index/round((length(song)-x(end))/spacingFactor))))
    n = x(end)+(index-1)*spacingFactor+1;
    % Generate random number
    if (rand < trainFrac)
        % Go to training set if below fracTrain
        trainAns(length(trainAns)+1) = song(n);
        for sub = 1:length(x)
            trainData(length(x)-sub+1,length(trainAns)) = song(n-x(length(x)-sub+1));
        end
    else
        % Else go to validation set
        validAns(length(validAns)+1) = song(n);
        for sub = 1:length(x)
            validData(length(x)-sub+1,length(validAns)) = song(n-x(sub));
        end
    end
%     figure(2)
%     subplot(2,2,1)
%     imagesc(trainData)
%     subplot(2,2,3)
%     imagesc(validData)
%     subplot(2,2,2)
%     plot(trainAns)
%     subplot(2,2,4)
%     plot(validAns)
%     drawnow
end

% Actually apply the data
for cycle = 1:trainingCycles
    % disp(w')
    
    % First check the error produced by the current weights on validData
    validErr = 0;
    for index = 1:length(validAns)
        validErr = validErr + (validAns(index)-(2./(1+exp(-1*(w*validData(:,index))))-1))^2;
    end
    
    errTrace(cycle) = validErr/length(validAns);
    % plot(errTrace),drawnow
    
    %See how bad the current weights do on the training set
    trainErr = 0;
    for index = 1:length(trainAns)
        trainErr = trainErr + (trainAns(index)-(2./(1+exp(-1*(w*trainData(:,index))))-1))^2;
    end
    
    % Record the direction each weight needs to move according to each
    % training example
    for layer = 1:layers
        for col = 1:length(weights{layer}(1,:))
            for row = 1:length(weights{layer}(:,1))
                %disp([cycle layer col row])
                % Check if increasing reduces error
                weights{layer}(row,col) = weights{layer}(row,col)+learnRate;
                wUP = weights{layers};
                for sub3 = 1:layers-1
                    wUP = wUP*weights{layers-sub3};
                end

                % Check if decreasing reduces error
                weights{layer}(row,col) = weights{layer}(row,col)-2*learnRate;
                wDN = weights{layers};
                for sub3 = 1:layers-1
                    wDN = wDN*weights{layers-sub3};
                end

                % Reset weights
                weights{layer}(row,col) = weights{layer}(row,col)+learnRate;
                
                errUP = 0;
                for index = 1:length(trainAns)
                    errUP = errUP + (trainAns(index)-(2./(1+exp(-1*(wUP*trainData(:,index))))-1))^2;
                end
                
                errDN = 0;
                for index = 1:length(trainAns)
                    errDN = errDN + (trainAns(index)-(2./(1+exp(-1*(wDN*trainData(:,index))))-1))^2;
                end

                if (min(errUP,errDN) < trainErr)% if new error smaller than err
                    if (errUP < errDN)
                        tempWeights{layer}(row,col) = tempWeights{layer}(row,col) + 1;
                    else % errDN < errUP
                        tempWeights{layer}(row,col) = tempWeights{layer}(row,col) - 1;
                    end
                % else neither of them improved the error, do nothing
                end     
            end
        end
    end
    % Finally, update the weights. *Hopefully* this leads to less error.
    for layer = 1:layers
        % Multiply tempWeights by the learning rate, and divide it by the
        % number of training examples to average out the impact of each
        weights{layer} = weights{layer} + tempWeights{layer}*learnRate/length(trainAns);
        % Reset tempWeights for next cycle
        tempWeights{layer} = tempWeights{layer}*0;
    end
    w = weights{layers};
    for index = 1:layers-1
        w = w*weights{layers-index};
    end
    disp(strcat('Training Progress:`',cycle/trainingCycles))
end

%% % Now compare the generated songs after training
prompt = nan(inputLength,1);
for index = 1:x(end)*3 %how many song datapoints to generate; Fs will make one second
    for subindex = 1:inputLength
        prompt(inputLength-subindex+1) = newSong(end - x(subindex));
    end
    newSong(length(newSong)+1) = 2./(1+exp(-1*(w*prompt)))-1;
end

figure
subplot(2,1,1)
plot(randSong)
subplot(2,1,2)
plot(newSong)

% sound((newSong(x(end):end)), Fs)