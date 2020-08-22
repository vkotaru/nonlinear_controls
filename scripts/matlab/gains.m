dt = 0.005;
A = [eye(3), dt*eye(3); zeros(3,3), eye(3)]
J = diag([112533, 36203, 42673])*1e-6;
B = [zeros(3,3); inv(J)]*dt
Q = eye(6); R = eye(3);
[K, P] = dlqr(A, B, Q, R)
Q = diag([1000, 1000, 1000, 100, 100, 100])
[K, P] = dlqr(A, B, Q, R)