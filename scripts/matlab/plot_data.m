dt = 1/200;
T = dt:dt:20;

close all;
figure;
hold on;
plot(T, afilempc.VarName17,'r');
plot(T, afilegeo.VarName21,'g');
latex_legend({'mpc','geo'});

figure;
hold on;
plot(T, afilempc.VarName10,'r');
plot(T, afilegeo.VarName10,'g');
latex_legend({'mpc','geo'});
latex_title('$$\eta_y$$');


figure;
hold on;
plot(T, afilempc.VarName13,'r');
plot(T, afilegeo.VarName13,'g');
latex_legend({'mpc','geo'});
latex_title('$$\Omega_y$$');



