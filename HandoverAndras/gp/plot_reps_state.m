function plot_reps_state(savem)
figure()
subplot(4,2,1), plot(savem(1:end, 1:4)), ylabel('mean')
subplot(4,2,2), plot(savem(2:end, 10)), ylabel('E[R]')
subplot(4,2,3), plot(savem(1:end, 5:6)), ylabel('std')
subplot(4,2,4), plot(savem(2:end, 7:9)), ylabel('eta,theta')
subplot(4,2,5), plot(savem(2:end, 12)), ylabel('Dkl')
subplot(4,2,6), plot(savem(1:end, 13)), ylabel('iDkl')
subplot(4,2,7), plot(savem(2:end, 11).^.5), ylabel('Var[R]^{0.5}')
