function NN_SimpleAdd()
    add = input('Additive factor = ');
    guessAdd = rand;

    count = 0;
    learn = 1;
    %learnRate = 1.05;
    error = 10^6;
    errVec = [];
    guessVec = [];
    figure
    while (error>1e-6)
        count = count + 1;

        inputData = rand(1,100);

        outputData = add+inputData;
        outputGuess = guessAdd+inputData;

        newAdd = randn*learn + guessAdd;
        newGuess = newAdd+inputData;

        currErr = sum((outputData - outputGuess).^2);
        newErr = sum((outputData - newGuess).^2);
        if (newErr < currErr)
            guessAdd = newAdd;
            error = newErr;
        else
            error = currErr;
        end

        errVec(count) = error;
        guessVec(count) = guessAdd;
        if (count>1000) && ((guessVec(count)-guessVec(count-1000))/guessVec(count)<0.001)
            break
        end
        
    end

    guessAdd
    count
    
    subplot(2,1,1)
    semilogx(guessVec)
    title('Guess')
    subplot(2,1,2)
    loglog(errVec)
    title('Error')
end