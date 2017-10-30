function NN_Simple_Mult()
    factor = input('Multiplicative factor = ');
    guessFactor = rand;

    count = 0;
    learn = 1;
    learnRate = 1.05;
    error = 10^6;
    errVec = [];
    while (error>1e-6)
        count = count + 1;

        inputData = rand(1,100);
        outputData = factor*inputData;

        outputGuess = guessFactor*inputData;

        newFactor = randn*learn + guessFactor;
        newGuess = newFactor*inputData;

        currErr = sum((outputData - outputGuess).^2);
        newErr = sum((outputData - newGuess).^2);
        if (newErr < currErr)
            guessFactor = newFactor;
            error = newErr;
        else
            error = currErr;
        end

        errVec(count) = error;
    end

    guessFactor
    count

    loglog(errVec)
end