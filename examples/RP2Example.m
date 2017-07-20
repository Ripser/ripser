%Programmer: Chris Tralie
%Purpose: To generate 600 random points on RP2 and show the difference
%between the H1/H2 diagrams with Z/2 and Z/3 coefficients 
addpath('..');
N = 200;
rng(N);
X = randn(N, 3);
X = bsxfun(@times, X, 1./sqrt(sum(X.^2, 2)));
D = acos(abs(X*X'));
D = real(D);
D = D + D';
D = D-diag(diag(D));

disp('Doing Z2...');
Is2 = ripserDM(D, 2, 2);
disp('Finished Z2');
disp('Doing Z3...');
Is3 = ripserDM(D, 3, 2);
disp('Finished Z3');

subplot(121);
plotDGM(Is2{2}, 'r', 20, 0), axis square;
hold on;
plotDGM(Is2{3}, [0, 0.5, 0]), axis square;
title('Z2 Coefficients');

subplot(122);
plotDGM(Is3{2}, 'r', 20, 0), axis square;
hold on;
plotDGM(Is3{3}, [0, 0.5, 0]),axis square;
title('Z3 Coefficients');