NN = 5;
lambda = zeros(NN,1);
perr = lambda;

for nn=1:NN
    BasisPerAE = nn;
    spectralFV;
    lambda(nn) = check;
    perr(nn) = errL2v;
end

plot(lambda/lambda(1),'*-');
hold on;
plot(perr/perr(1),'ro--');
hold off;