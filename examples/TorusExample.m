%Programmer: Jose Perea
%Purpose: To generate 400 random points on the 2-D torus and show the difference
%between the H1/H2 diagrams with Z/2 and Z/3 coefficients 
addpath('..');
N = 400; 
rng(N);
angle = 2*pi*rand(2,N);

data = [cos(angle(1,:)); sin(angle(1,:)); cos(angle(2,:)); sin(angle(2,:))]';

% geodesic
Dgeo = sqrt(acos(data(:,1:2)*data(:,1:2)').^2 + acos(data(:,3:4)*data(:,3:4)').^2);
Dgeo = Dgeo - diag(diag(Dgeo));

% ambient
%Damb = squareform(pdist(data));
XSqr = sum(data.^2, 2);
Damb = bsxfun(@plus, XSqr(:), XSqr(:)') - 2*(data*data');
Damb = 0.5*(Damb + Damb');
Damb(Damb < 0) = 0;

disp('Doing Z2...');
Is2geo = ripserDM(Dgeo, 2, 2);
Is2amb = ripserDM(Damb, 2, 2);
disp('Finished Z2');
disp('Doing Z3...');
Is3geo = ripserDM(Dgeo, 3, 2);
Is3amb = ripserDM(Damb, 3, 2);
disp('Finished Z3');

figure
subplot(221);
plotDGM(Is2geo{2}, 'r', 20, 0), axis square;
hold on;
plotDGM(Is2geo{3}, [0, 0.5, 0]), axis square;
title('Z2 Coefficients / Geodesic');

subplot(222);
plotDGM(Is3geo{2}, 'r', 20, 0), axis square;
hold on;
plotDGM(Is3geo{3}, [0, 0.5, 0]),axis square;
title('Z3 Coefficients / Geodeisc');

subplot(223);
plotDGM(Is2amb{2}, 'r', 20, 0), axis square;
hold on;
plotDGM(Is2amb{3}, [0, 0.5, 0]), axis square;
title('Z2 Coefficients / Ambient');

subplot(224);
plotDGM(Is3amb{2}, 'r', 20, 0), axis square;
hold on;
plotDGM(Is3amb{3}, [0, 0.5, 0]),axis square;
title('Z3 Coefficients / Ambient');