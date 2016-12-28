clear all, close all,

addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

ridge = 1e-4;

fixedActivation = 0.2;

userData;
figure
for j =1:length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    
   ixAbs = data.absFeedback(:, 1);
   normSamples = bsxfun(@rdivide, bsxfun(@minus, data.samples(ixAbs, :), data.polparMinLimit), data.polparMaxLimit-data.polparMinLimit);
   absFeedback = data.absFeedback(:, 2);
   
   normLearned = bsxfun(@rdivide, bsxfun(@minus, data.policyMean(end, :), data.polparMinLimit), data.polparMaxLimit-data.polparMinLimit);
   normInit = bsxfun(@rdivide, bsxfun(@minus, data.policyMean(1, :), data.polparMinLimit), data.polparMaxLimit-data.polparMinLimit);
   diffLearned = abs(maha(normLearned, normInit)).^.5;
   
   rho = [];
   for i = 10:length(ixAbs)
       
       dist = abs(maha(normSamples(1:i, :), normSamples(1:i, :))).^.5;
       dist1 = abs(maha(normSamples(1:i, 1), normSamples(1:i, 1))).^.5;
       dist2 = abs(maha(normSamples(1:i, 2), normSamples(1:i, 2))).^.5;
       dist3 = abs(maha(normSamples(1:i, 3), normSamples(1:i, 3))).^.5;
       dist4 = abs(maha(normSamples(1:i, 4), normSamples(1:i, 4))).^.5;
       
       rewDist = abs(maha(absFeedback(1:i), absFeedback(1:i))).^.5;
       dummy = [];
       dummy2 = [];
       dummy11 = [];
       dummy12 = [];
       dummy13 = [];
       dummy14 = [];
       
       
       for k = 1:i
           
           dummy(k, :) = dist(k, [1:(k-1), (k+1):i]);
           dummy11(k, :) = dist1(k, [1:(k-1), (k+1):i]);
           dummy12(k, :) = dist2(k, [1:(k-1), (k+1):i]);
           dummy13(k, :) = dist3(k, [1:(k-1), (k+1):i]);
           dummy14(k, :) = dist4(k, [1:(k-1), (k+1):i]);
           
           
           dummy2(k, :) = rewDist(k, [1:(k-1), (k+1):i]);
           dummy3 = corrcoef(dummy(k, :)', dummy2(k, :)');
             rho(i, k) = dummy3(1, 2);
             
             dummy31 = corrcoef(dummy11(k, :)', dummy2(k, :)');
             dummy32 = corrcoef(dummy12(k, :)', dummy2(k, :)');
             dummy33 = corrcoef(dummy13(k, :)', dummy2(k, :)');
             dummy34 = corrcoef(dummy14(k, :)', dummy2(k, :)');
             rho(i, k) = dummy3(1, 2);
             rho1(i, k) = dummy31(1, 2);
             rho2(i, k) = dummy32(1, 2);
             rho3(i, k) = dummy33(1, 2);
             rho4(i, k) = dummy34(1, 2);
             
       end
       meanRho(i) = mean(rho(i, 1:i));
       stdRho(i) = std(rho(i, 1:i));
       
       meanRho1(i) = mean(rho1(i, 1:i));
       stdRho1(i) = std(rho1(i, 1:i));
       
       meanRho2(i) = mean(rho2(i, 1:i));
       stdRho2(i) = std(rho2(i, 1:i));
       
       
       meanRho3(i) = mean(rho3(i, 1:i));
       stdRho3(i) = std(rho3(i, 1:i));
       
       meanRho4(i) = mean(rho4(i, 1:i));
       stdRho4(i) = std(rho4(i, 1:i));
       
       
       dist = dummy;
       rewDist = dummy2;
       
   end
   
   clf, subplot(2,2,1), plot(meanRho), hold on, plot(stdRho, '--')
   title([user.names{j}, ', mean(meanRho) = ', num2str(mean(meanRho))]);
   subplot(2,2,2), plot(rho(end, :))
   title([user.names{j}, ', mean(corrcoef All) = ', num2str(mean(rho(end, :))), ' diffLearned = ', num2str(diffLearned)]);
   subplot(2,2,3), plot([meanRho1; meanRho2; meanRho3; meanRho4]'), ylabel('meanRho Pars')
   subplot(2,2,4), plot([stdRho1; stdRho2; stdRho3; stdRho4]'), ylabel('meanRho Pars')
   drawnow
   
   pause

end