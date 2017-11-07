%@@TODO: Output weights over coarse digitization instead of single number?
%@@TODO: Check if normalizing the song actually helps meaningfully.
%@@TODO: Test batch vs online updating. Note that online is going to
%        require randomly picking from among training examples to avoid
%        biasing the results by always feeding in early examples first.

clear all
close all

% Set parameters
compressionFactor = 2; % Retain 1/cF of the song data
secondsBack = 10; % Length of time to be used as input
learnRate = 0.001; % Smaller = slower adjustment of weights
trainFrac = 0.8; % Portion of data to be devoted to training set
trainingCycles = 100; % # of times training data is run through 
    % Currently arbitrary - find optimal number by looking at error trace
errTrace = zeros(1,trainingCycles); % Store error on validataion set
spacingFactor = 100; % Space between datapoints used for training

% First grab the song. Location will vary for your computer, obvs.
% I *think* audioread assumes an nx2 mp3 format...?
[song,Fs] = audioread('C:\Users\Oli\Downloads\Simple Gifts.mp3');

% Make it one-channel
song = mean(song,2);

% 'Simple Gitfts' has a chunk in front of zeros, discard them
song = song(576:end);

% Lossy compression by picking out every (compressionFactor)-th data point
song = song(1:compressionFactor:end);
Fs = Fs/compressionFactor;

% sound(song,Fs)
% clear sound
disp('Song Obtained')
%%
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
%%
%Random initial weights on (-valueMax,valueMax)
valueMax = 1;
weights{1} = 2*valueMax*rand(round(inputLength/2),inputLength)-valueMax;
biases{1} = 2*valueMax*rand(ceil(length(x)/2),1)-valueMax;
tempWeights{1} = weights{1}*0;
tempBiases{1} = biases{1}*0;
for index = 2:layers
    weights{index} = 2*valueMax*(rand(ceil(length(weights{index-1}(:,1))/2),ceil(length(weights{index-1}(1,:))/2)))-valueMax;
    tempWeights{index} = weights{index}*0;
    
    biases{index} = 2*valueMax*(rand(ceil(length(biases{index-1})/2),1))-valueMax;
    tempBiases{index} = biases{index}*0;
end
disp('Initialized Weights and Biases')

%%

% Create an example song w/ the random weights
% Choose random chunk of song as a seed.

close all
disp('Producing Random Song')
numNewSounds = x(end); % How many audi datapoints to generate
m = randi([x(end) + 1, length(song)]);
newSong = repmat([song(m-x(end) : m) ; nan(numNewSounds,1)],1,trainingCycles+1);

for index = 1:numNewSounds
    newPoint = nan(inputLength,1);
    for subindex = 1:inputLength
        newPoint(inputLength-subindex+1) = newSong(x(end)+1+index - x(subindex));
    end
    
    for subindex = 1:layers
        newPoint = weights{subindex}*newPoint + biases{subindex};
        newPoint = 2./(1+exp(-1*newPoint))-1;
    end
    
    newSong(x(end)+1+index,1) = newPoint;
end
figure
imagesc(newSong)
figure
plot(newSong(:,1))

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
    for example = 1:length(validAns)
        newPoint = validData(:,example);
        for subindex = 1:layers
           newPoint = weights{subindex}*newPoint + biases{subindex};
           newPoint = 2./(1+exp(-1*newPoint))-1;
        end
        validErr = validErr + (validAns(example)-newPoint)^2;
    end
    errTrace(cycle) = validErr/length(validAns);
    figure(10)
    plot(errTrace)
    title('Mean Validation Error'),drawnow
    
    % Record the direction each weight and bias needs to move according to
    % each training example
    for example = 1:length(trainAns)
        if mod(example,10) == 0
            disp(strcat(strcat('Training:`',num2str(cycle - 1 + example/length(trainAns))),strcat('/',num2str(trainingCycles))))
        end
        % check original error
        newPoint = trainData(:,example);
        for subindex = 1:layers
           newPoint = weights{subindex}*newPoint + biases{subindex};
           newPoint = 2./(1+exp(-1*newPoint))-1;
        end
        err0 = (trainAns(example)-newPoint)^2;
        
        for layer = 1:layers
            for n = 1:length(biases{layer})
                biases{layer}(n) = biases{layer}(n)+learnRate;
                % Check error moving up and down in a small step
                newPoint = trainData(:,example);
                for subindex = 1:layers
                   newPoint = weights{subindex}*newPoint + biases{subindex};
                   newPoint = 2./(1+exp(-1*newPoint))-1;
                end
                errU = (trainAns(example)-newPoint)^2;
                
                biases{layer}(n) = biases{layer}(n)-2*learnRate;
                newPoint = trainData(:,example);
                for subindex = 1:layers
                   newPoint = weights{subindex}*newPoint + biases{subindex};
                   newPoint = 2./(1+exp(-1*newPoint))-1;
                end
                errD = (trainAns(example)-newPoint)^2;
                
                [temp,location] = min([err0,errU,errD]);
                switch location
                case 1
                    break
                case 2
                    tempBiases{layer}(n) = tempBiases{layer}(n)+1;
                case 3
                    tempBiases{layer}(n) = tempBiases{layer}(n)-1;
                otherwise
                    disp('Something is very wrong!')
                    pause(10000)
                    break
                end
                % Reset the bias (+ -- + => 0)
                biases{layer}(n) = biases{layer}(n)+learnRate;
            end
            for col = 1:size(weights{layer},2)
                for row = 1:size(weights{layer},1)
                    weights{layer}(row,col) = weights{layer}(row,col)+learnRate;
                    newPoint = trainData(:,example);
                    for subindex = 1:layers
                       newPoint = weights{subindex}*newPoint + biases{subindex};
                       newPoint = 2./(1+exp(-1*newPoint))-1;
                    end
                    errU = (trainAns(example)-newPoint)^2;

                    weights{layer}(row,col) = weights{layer}(row,col)-2*learnRate;
                    newPoint = trainData(:,example);
                    for subindex = 1:layers
                       newPoint = weights{subindex}*newPoint + biases{subindex};
                       newPoint = 2./(1+exp(-1*newPoint))-1;
                    end
                    errD = (trainAns(example)-newPoint)^2;

                    % Reset the weight                
                    weights{layer}(row,col) = weights{layer}(row,col)+learnRate;

                    [temp,location] = min([err0,errU,errD]);
                    switch location
                        case 1
                            break
                        case 2
                            tempWeights{layer}(row,col) = tempWeights{layer}(row,col)+1;
                        case 3
                            tempWeights{layer}(row,col) = tempWeights{layer}(row,col)-1;
                        otherwise
                            disp('Something is very wrong!')
                            pause(10000)
                            break
                    end
                end
            end  
        end
    end
    
    % Update the weights and biases
    % *Hopefully* this leads to less error.
    for layer = 1:layers
        % Multiply tempWeights by the learning rate, and divide it by the
        % number of training examples to average out the impact of each
        weights{layer} = weights{layer} + tempWeights{layer}*learnRate/length(trainAns);
        % Reset tempWeights for next cycle
        tempWeights{layer} = tempWeights{layer}*0;
        
        biases{layer} = biases{layer} + tempBiases{layer}*learnRate/length(trainAns);
        tempBiases{layer} = tempBiases{layer}*0;
    end
    
    % Create a new song based on the new weights and biases, with the same
    % seed as the initial random song.
    for index = 1:numNewSounds
        newPoint = nan(inputLength,1);
        for subindex = 1:inputLength
            newPoint(inputLength-subindex+1) = newSong(x(end)+1+index - x(subindex));
        end

        for subindex = 1:layers
            newPoint = weights{subindex}*newPoint + biases{subindex};
            newPoint = 2./(1+exp(-1*newPoint))-1;
        end

        newSong(x(end)+1+index,cycle+1) = newPoint;
    end
end

%% % Now compare the generated songs after training

figure
imagesc(newSong)
figure
for plt = 1:trainingCycles+1
    subplot(1,trainingCycles+1,plt)
    plot(newSong(:,plt))
end