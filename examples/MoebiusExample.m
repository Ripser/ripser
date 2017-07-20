%Moebius strip example.  A fake H1 class exists with Z2 coefficients
%which does not exist with Z3 coefficients
addpath('..');

t = linspace(0, 2*pi, 500);
t = t(:);
X = [(2 + cos(t)).*cos(2*t) (2 + cos(t)).*sin(2*t) sin(t)];


disp('Doing Z2...');
Is2 = ripserPC(X, 2, 1);
disp('Finished Z2');
disp('Doing Z3...');
Is3 = ripserPC(X, 3, 1);
disp('Finished Z3');

subplot(121);
plotDGM(Is2{2}, 'r', 20);
axis square;
title('Z2 Coefficients');

subplot(122);
plotDGM(Is3{2}, 'r', 20);
axis square;
title('Z3 Coefficients');
