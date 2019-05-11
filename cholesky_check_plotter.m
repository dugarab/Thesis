c_before = ((1/rc)^0.5)*ones(1,total_iters);
c_after = ((1/rc)*condition_num_gain).^0.5;

gain_calculated = (log(c_after+1) - log(c_after-1))./(log(c_before+1) - log(c_before-1));
t = 1:total_iters;
figure(1);
plot(t, gain_calculated, '.',t, iteration_gain_cg,'x',t,iteration_gain_gmres,'d');
legend("Theoretical gain","observed CG gain", "observed GMRES gain");
ylabel("Gain");
xlabel('Sample number');
title(['Comparison of gain from cholesky N = ' num2str(500) ', cond = ' num2str(1/rc)]);

figure(2);

semilogy(t, residue_cg1,'.', t, residue_cg2,'x');
legend("x^TAx rel Error","b-Ax rel Error");
ylabel("Error");
xlabel('Sample number');
title("Comparison of errors from conjugate gradient");

figure(3);
plot(t, total_iters_cg,'.', t, total_iters_gmres,'x');
legend("iters by CG","iters by GMRES");
ylabel("Iteration Count");
xlabel('Sample number');
title("Comparison of non-preconditioned iteration count by CG and GMRES");
