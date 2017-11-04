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
spacingFactor = 1000;

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

song = 2*(song-minSong)/(maxSong-minSong)-1;

% sound(song,Fs)
% clear sound
disp('Song Obtained')
%%
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

weights{1} = 2*rand(round(inputLength/2),inputLength)-1;
biases{1} = 2*rand(ceil(length(x)/2),1)-1;
tempWeights{1} = weights{1}*0;
tempBiases{1} = biases{1}*0;
for index = 2:layers
    weights{index} = 2*(rand(ceil(length(weights{index-1}(:,1))/2),ceil(length(weights{index-1}(1,:))/2)))-1;
    tempWeights{index} = weights{index}*0;
    
    biases{index} = 2*(rand(ceil(length(biases{index-1})/2),1))-1;
    tempBiases{index} = biases{index}*0;
end
disp('Initialized Weights and Biases')
%%
% Create an example song w/ the random weights
% Choose random chunk of song as a seed.
disp('Producing Random Song')
m = randi([x(end) + 1, length(song)]);
randSong = song(m-x(end) : m);
 
newSong = randSong; %This will be used to compare outputs later
for index = 1:x(end)*3 %how many song datapoints to generate
    newPoint = nan(inputLength,1);
    for subindex = 1:inputLength
        newPoint(inputLength-subindex+1) = randSong(end - x(subindex));
    end
    
    for subindex = 1:layers
        newPoint = weights{subindex}*newPoint + biases{subindex};
        newPoint = 2./(1+exp(-1*newPoint))-1;
    end
    
    randSong(length(randSong)+1) = newPoint;
end

%%
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

%%
% Actually apply the data
for cycle = 1:trainingCycles    
    % First check the error produced by the current weights on validData
    validErr = 0;
    for index = 1:length(validAns)
        newPoint = nan(inputLength,1);
        for subindex = 1:inputLength
            newPoint(inputLength-subindex+1) = validData(end - subindex + 1,index);
        end
        for subindex = 1:layers
           newPoint = weights{subindex}*newPoint + biases{subindex};
           newPoint = 2./(1+exp(-1*newPoint))-1;
        end
        validErr = validErr + (validAns(index)-newPoint)^2;
    end    
    errTrace(cycle) = validErr/length(validAns);
    plot(errTrace)
    title('Validation Error'),drawnow
    
    %See how bad the current weights do on the training set
    trainErr = 0;
    for index = 1:length(trainAns)
        newPoint = nan(inputLength,1);
        for subindex = 1:inputLength
            newPoint(inputLength-subindex+1) = trainData(end - subindex + 1);
        end
        for subindex = 1:layers
           newPoint = weights{subindex}*newPoint + biases{subindex};
           newPoint = 2./(1+exp(-1*newPoint))-1;
        end
        trainErr = trainErr + (trainAns(index)-newPoint)^2;
    end
    
    % Record the direction each weight and bias needs to move according to
    % each training example
    for layer = 1:layers
        for n = 1:length(biases{layer})
            biases{layer}(n) = biases{layer}(n)+learnRate;
            % Re-do the whole error calculation
            trainErrTemp = 0;
            for index = 1:length(trainAns)
                newPoint = nan(inputLength,1);
                for subindex = 1:inputLength
                    newPoint(inputLength-subindex+1) = trainData(end - subindex + 1);
                end
                for subindex = 1:layers
                   newPoint = weights{subindex}*newPoint + biases{subindex};
                   newPoint = 2./(1+exp(-1*newPoint))-1;
                end
                trainErrTemp = trainErrTemp + (trainAns(index)-newPoint)^2;
            end
            if trainErrTemp < trainErr
                tempBiases{layer}(n) = tempBiases{layer}(n)+1;
            end
            % Repeat for reducing that bias
            biases{layer}(n) = biases{layer}(n)-2*learnRate;
            % Re-do the whole error calculation
            trainErrTemp = 0;
            for index = 1:length(trainAns)
                newPoint = nan(inputLength,1);
                for subindex = 1:inputLength
                    newPoint(inputLength-subindex+1) = trainData(end - subindex + 1);
                end
                for subindex = 1:layers
                   newPoint = weights{subindex}*newPoint + biases{subindex};
                   newPoint = 2./(1+exp(-1*newPoint))-1;
                end
                trainErrTemp = trainErrTemp + (trainAns(index)-newPoint)^2;
            end
            if trainErrTemp < trainErr
                tempBiases{layer}(n) = tempBiases{layer}(n)-1;
            end
            % Reset the bias; + -- + => 0
            biases{layer}(n) = biases{layer}(n)+learnRate;
        end
        for col = 1:length(weights{layer}(1,:))
            for row = 1:length(weights{layer}(:,1))
                % Check if increasing reduces error
                weights{layer}(row,col) = weights{layer}(row,col)+learnRate;
                % Re-do the whole error calculation
                trainErrTemp = 0;
                for index = 1:length(trainAns)
                    newPoint = nan(inputLength,1);
                    for subindex = 1:inputLength
                        newPoint(inputLength-subindex+1) = trainData(end - subindex + 1);
                    end
                    for subindex = 1:layers
                       newPoint = weights{subindex}*newPoint + biases{subindex};
                       newPoint = 2./(1+exp(-1*newPoint))-1;
                    end
                    trainErrTemp = trainErrTemp + (trainAns(index)-newPoint)^2;
                end
                if trainErrTemp < trainErr
                    tempWeights{layer}(row,col) = tempWeights{layer}(row,col)+1;
                end
                % Check if reducing reduces error
                weights{layer}(row,col) = weights{layer}(row,col)-2*learnRate;
                % Re-do the whole error calculation
                trainErrTemp = 0;
                for index = 1:length(trainAns)
                    newPoint = nan(inputLength,1);
                    for subindex = 1:inputLength
                        newPoint(inputLength-subindex+1) = trainData(end - subindex + 1);
                    end
                    for subindex = 1:layers
                       newPoint = weights{subindex}*newPoint + biases{subindex};
                       newPoint = 2./(1+exp(-1*newPoint))-1;
                    end
                    trainErrTemp = trainErrTemp + (trainAns(index)-newPoint)^2;
                end
                if trainErrTemp < trainErr
                    tempWeights{layer}(row,col) = tempWeights{layer}(row,col)-1;
                end
                % Reset the weight
                weights{layer}(row,col) = weights{layer}(row,col)+learnRate;
                end     
            end
    end

    % Finally, update the weights and biases
    % *Hopefully* this leads to less error.
    for layer = 1:layers
        % Multiply tempWeights by the learning rate, and divide it by the
        % number of training examples to average out the impact of each
        weights{layer} = weights{layer} + tempWeights{layer}*learnRate;%/length(trainAns);
        % Reset tempWeights for next cycle
        tempWeights{layer} = tempWeights{layer}*0;
        
        biases{layer} = biases{layer} + tempBiases{layer}*learnRate;%/length(trainAns);
        tempBiases{layer} = tempBiases{layer}*0;
    end
    disp(strcat('Training Progress:`',num2str(cycle/trainingCycles)))
end

%% % Now compare the generated songs after training
newPoint = nan(inputLength,1);
for index = 1:x(end)*3 %how many song datapoints to generate
    newPoint = nan(inputLength,1);
    for subindex = 1:inputLength
        newPoint(inputLength-subindex+1) = newSong(end - x(subindex));
    end
    
    for subindex = 1:layers
        newPoint = weights{subindex}*newPoint + biases{subindex};
        newPoint = 2./(1+exp(-1*newPoint))-1;
    end
    
    newSong(length(newSong)+1) = newPoint;
end

figure
subplot(2,1,1)
plot(randSong)
subplot(2,1,2)
plot(newSong)

% sound((newSong(x(end):end)), Fs)