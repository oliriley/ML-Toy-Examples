function NN_SimpleMult()
    factor = input('Multiplicative factor = ');
    guessFactor = rand;

    count = 0;
    learn = 1;
    learnRate = 1.05;
    error = 10^6;
    errVec = [];
    guessVec = [];
    figure
    while (error>1e-6)
        count = count + 1;

        inputData = rand(1,100);
        outputData = factor*inputData;

        outputGuess = guessFactor*inputData;

        newFactor = randn*learn + guessFactor;
        newGuess = newFactor*inputData;

        currErr = sum((outputData - outputGuess).^2);
        newErr = sum((outputData - newGuess).^2);
        if (newErr <= currErr)
            guessFactor = newFactor;
            error = newErr;
        else
            error = currErr;
        end

        errVec(count) = error;
        guessVec(count) = guessFactor;
        if (count>1000) && ((guessVec(count)-guessVec(count-1000))/guessVec(count)<0.001)
            break
        end
    end

    guessFactor
    count
    
    subplot(2,1,1)
    semilogx(guessVec)
    title('Multiplier Guess')
    subplot(2,1,2)
    loglog(errVec)
    title('Squared Error')
end