% hard example of using rmvnrnd
S = [ 1 0.98; 0.98 1 ];
m = [ -3; 0 ];

% draw a sample of 10000 from the unrestricted MVN
X1 = rmvnrnd(m, S, 10000);

% draw a sample of 1000 from the unit square
A = [-eye(2); eye(2)];
b = [0;0;1;1];
[X2, rho2, nar2, ngibbs2] = rmvnrnd(m, S, 1000, A, b);

% draw a sample of 2000 from a half plane
A = [ 1 -1];
b = [ -4];
[X3, rho3, nar3, ngibbs3] = rmvnrnd(m, S, 2000, A, b);

% Plot the points
clf
plot(X1(:,1),X1(:,2),'g.');
hold on
grid on
plot(X2(:,1),X2(:,2),'b.');
plot(X3(:,1),X3(:,2),'r.');

axis equal
legend('unrestricted','unit square','x - y <= -4','Location','SouthEast')