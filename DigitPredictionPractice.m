%DigitPredictionPractice
fid0 = fopen('data0.txt','r');
fid1 = fopen('data1.txt','r');

data0 = zeros(28,28,1000);
data1 = zeros(28,28,1000);

for index = 1:1000
    data0(:,:,index) = round((fread(fid0,[28 28],'uchar')+1)/256);
    data1(:,:,index) = round((fread(fid1,[28 28],'uchar')+1)/256);
end

%Derive input features
inputs = zeros(2000,3);
