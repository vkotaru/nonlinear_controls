N = 5;
nx = 6;
nu = 3;
dt = 1/200;

A = eye(nx);
A(1:3, 4:end) = A(1:3, 4:end) + dt*eye(3);
B = [0.5*dt*dt*eye(3); dt*eye(3)];
Q = 0.1*eye(6);
R = 0.01*eye(3);
P = 0.1*eye(6);


Sx = zeros(nx*(N+1), nx);
Su = zeros(nx*(N+1), N*nu);
Sx(1:nx, :) = eye(nx);
Qbar = [];

for i = 1:N
    Sx(nx*(i)+1:nx*(i+1),:) = Sx(nx*(i-1)+1:nx*(i),:)*A;
    if (i==1)
        tmp = B;
    else
       tmp = A*Su(nx*(i-1)+1:nx*(i),1:nu);
    end
    Su(nx*(i)+1:nx*(i+1),:) = [tmp, Su(nx*(i-1)+1:nx*(i),1:end-nu)];
    
    Qbar= blkdiag(Qbar, Q);
end
Qbar= blkdiag(Qbar, P);


