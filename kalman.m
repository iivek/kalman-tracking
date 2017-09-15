%% A small example of tracking using Kalman filtering

% Let's create an artificial trajectory in 2D
T = pi/16;
t = [pi/16:T:8*pi];
x = -cos(t/4);
y = sin(t);
% add observation noise
x = x + normrnd(0,0.05,size(x));
y = y + normrnd(0,0.01,size(x));


hold on, plot(x,y, 'b')
plot(x,y, '.')
axis(1.2*[-1 1 -1 1])


% Let's assume a model of movement where speed and acceleration are taken
% into account and expressed as follows:
% v[k] = (x[k]-x[k-1])/T; % first discrete derivative - speed
% a[k] = (v[k]-v[k-1])/T = (x([k]-2*x[k-1]+x[k-2])/T^2;  % second discrete derivative - acceleration
% x[k+1] = x[k] + T*v[k] + [T^2]/2 * a[k]
% For the y-coordinate, the expressions are the same.

% Transition matrix A, based on the above expressions
A = [5/2 -2 1/2; 1 0 0; 0 1 0];  
A = [A, zeros(size(A)); zeros(size(A)), A];

% The only available measurements are position-x and position-y
H = zeros(2, size(A,2));
H(1,1) = 1;
H(2,size(A,2)/2+1) = 1;

% Kalman filtering initialization
P = eye(size(A));
sigma_W = 0.01; % process noise cov, let's say it's isotropic
Rww = sigma_W*eye(size(A));
Rvv = [0.05 0;
       0 0.01]; % observation noise covariance
X = zeros(size(A,1),1);
X(1) = x(1);
X(size(A,1)/2+1) = y(1);
Predikcija = zeros(2, length(t));

% Filtering
for i=2:length(t)
    % Prediction equations
    X = A*X;
    P = A*P*A' + Rww;
    % Correction equations
    Z = [x(i); y(i)];
    K = P*H' * inv(H*P*H' + Rvv);
    X = X + K*(Z-H*X);
    P = (eye(size(A))-K*H)*P;
    %
    prediction(:,i) = H*X;
end

plot(prediction(1,2:end), prediction(2,2:end), 'r')
plot(prediction(1,2:end), prediction(2,2:end), '.r')
hold off