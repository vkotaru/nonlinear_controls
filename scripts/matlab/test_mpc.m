N = 5;
nx = 6;
nu = 3;
dt = 1/200;

A = eye(nx);
A(1:3, 4:end) = A(1:3, 4:end) + dt*eye(3);
B = [0.5*dt*dt*eye(3); dt*eye(3)];

Q = blkdiag(100*eye(3), 20*eye(3));
R = eye(3);
[K, P, ~] = dlqr(A, B, Q, R)

% 
% Sx = zeros(nx*(N+1), nx);
% Su = zeros(nx*(N+1), N*nu);
% Sx(1:nx, :) = eye(nx);
% Qbar = [];
% Rbar = [];
% 
% xlb = 
% 
% for i = 1:N
%     Sx(nx*(i)+1:nx*(i+1),:) = Sx(nx*(i-1)+1:nx*(i),:)*A;
%     if (i==1)
%         tmp = B;
%     else
%        tmp = A*Su(nx*(i-1)+1:nx*(i),1:nu);
%     end
%     Su(nx*(i)+1:nx*(i+1),:) = [tmp, Su(nx*(i-1)+1:nx*(i),1:end-nu)];
%     
%     Qbar= blkdiag(Qbar, Q);
%     Rbar = blkdiag(Rbar, R);
% end
% Qbar= blkdiag(Qbar, P);
% 
% 
% lbA = 
% 
% 
