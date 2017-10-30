function NN_AddMult()
    mult = input('Multiplicative factor = ');
    add = input('Multiplicative factor = ');
    guessAdd = rand;
    guessMult = rand;

    count = 0;
    learn = 1;
    error = 10^6;
    errVec = [];
    addVec = [];
    multVec = [];
    figure
    while (error>1e-6)
        count = count + 1;

        inputData = rand(1,100);
        
        outputData = mult*inputData + add;
        outputGuess = guessMult*inputData + guessAdd;

        newMult = randn*learn + guessMult;
        newGuess = newMult*inputData + guessAdd;

        currErr = sum((outputData - outputGuess).^2);
        newErr = sum((outputData - newGuess).^2);
        if (newErr < currErr)
            guessMult = newMult;
        end
        
        newAdd = randn*learn + guessAdd;
        newGuess = guessMult*inputData + newAdd;

        currErr = sum((outputData - outputGuess).^2);
        newErr = sum((outputData - newGuess).^2);
        if (newErr < currErr)
            guessAdd = newAdd;
            error = newErr;
        else
            error = currErr;
        end

        errVec(count) = error;
        addVec(count) = guessAdd;
        multVec(count) = guessMult;
        
        if (count>1000) && ((multVec(count)-multVec(count-1000))/multVec(count)<0.001) && ((addVec(count)- addVec(count-1000))/ addVec(count)<0.001)
            break
        end
    end

    guessMult
    guessAdd
    count

    subplot(2,1,1)
    semilogx(multVec),hold on
    semilogx(addVec,'r')
    legend('Multiplier','Addition')
    title('Guesses')
    subplot(2,1,2)
    loglog(errVec)
    title('Squared Error')
end