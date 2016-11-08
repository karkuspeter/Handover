clear all, close all

addpath('../gp')
addpath('../gpml/')
startup
addpath('..')

[fPrior, prefs, context, samples, forceJerk, ixSlow] = humanFeedbackRobotHandover();
% sampleMean = mean(samples, 1);
% sampleStd = std(samples);
% samples = bsxfun(@minus, samples, mean(samples, 1));
% samples = bsxfun(@rdivide, samples, std(samples));
% context_orig = context;
% contextMean = mean(context);
% contextStd = std(context);
% context = context - mean(context);
% context =  context./std(context);
% fPrior(:, 2) = fPrior(:, 2) - 5;
% fPriorMean = mean(fPrior(:, 2));

% ixOkPrior = fPrior(:, 1) > 85;
% ixOkPref = and(prefs(:, 1) > 85, prefs(:, 2) > 85);
% 
% fPrior = fPrior(ixOkPrior, :);
% fPrior(:, 1) = fPrior(:, 1) - 85;
% prefs = prefs(ixOkPref, :);
% prefs = prefs -85;
% samples = samples(86:end, :);
% context = context(86:end);



sigp = 2;
sigf = 3;
siga = 1;
ridge = 1e-4;
% 
% loghyp = [log(std([samples, context])/2), log([sigf sigp ])]';
% loghyp_init = loghyp;
% 
%  opts.MaxFunEvals = 200;
%  opts.SaveVariables = 'off';
%  opts.TolFun = .1;
% 
%  foptsave = [];
%  lhoptsave = [];
% 
%  objfun = @(lh) pref_loghyp_opt(lh, [samples, context], prefs, fPrior, ridge, siga);
%  parfor i = 1:8
%      loghyp = [log(std([samples, context])), log([sigf sigp ])]';
%      [~, ~, ~, ~, out] = cmaes('pref_loghyp_opt', loghyp, ones(length(loghyp), 1), opts, [samples, context], prefs, fPrior, ridge, siga);
% 
% 
%      foptsave(i) = out.solutions.bestever.f;
%      lhoptsave(i, :) = out.solutions.bestever.x';
% 
%  end
%  foptsave(end+1) = objfun(loghyp_init);
%  lhoptsave( end+1, :) = loghyp_init(:)';
%  [foptsave; lhoptsave']
%  keyboard
% 
% ferike = [   121.1964  122.1184  121.8093  114.5125  123.4895  123.5342  121.5580  121.7061  120.1824  123.3770  136.8106
%     0.1849   -3.0063   -0.4322   -0.2636   -3.4790    0.6087    1.6781    2.7074    0.2059    4.0194   -0.6931
%     0.6246    0.4800    1.9497    1.8383   -0.3493    0.7296    0.4893    0.3006    0.7131   -2.8376   -0.6931
%     -0.4888    0.2223    0.1629    1.8550   -1.3833   -3.5704   -0.0077   -0.0948   -0.3338    0.5730   -0.6931
%     -0.5607    1.2466   -0.2528    0.4873    0.6609   -1.9853    1.0693    0.1824    0.9075   -0.1886   -0.6931
%     2.1370    1.1075    1.1803   -2.5439    1.2558   -5.0436   -2.6799   -0.0422   -0.1737    3.4791   -0.6931
%     1.3173   -0.8860    0.0740    4.8230    0.3680   -4.0459    1.0755   -0.9371   -0.0762   -0.0874   -0.6931
%     0.0292    2.6149   -2.8210    2.3531   -0.6038   -0.9538    2.9273    0.4117    2.6046    2.7481   -0.6931
%     -0.3690   -1.4591    2.0172    1.2893    1.0021   -5.0440   -0.8565   -1.8701   -2.2455   -2.1382   -0.6931
%     0.8829    0.6490    0.5641    0.8806    0.6280    0.6242    0.5744    0.6945    0.5962    0.5880    1.0986
%     -1.6610   -1.0627   -1.6408   -1.6353   -1.1450   -1.1101   -0.3616   -1.2880   -1.4350   -1.4074    0.6931];

% ferikeNoScale = [156.5111  141.0883  143.2487  142.1414  155.2194  147.1008  152.8237  144.3995  148.3330  144.4190  259.8900
%     1.2588    0.9325    0.4203    0.9354    3.7017    1.1848    1.8162    0.8092    0.7805    1.1569   -0.4938
%     5.8514    5.9334    6.2784    6.9072    5.9046    6.9972    6.2621    6.4987    5.3718    8.0854    4.2097
%     6.7374    6.7413    6.8178    6.3018    6.8772    6.1886    5.7832    6.2013    5.4376    7.1054    4.5829
%     5.1979    6.2801    5.3584    6.6772    6.9057    5.0858    4.4546    6.7051    6.8671    5.9147    4.1260
%     2.0895   -0.1099   -0.3390   -0.0790    2.0245    1.9401    4.0078    0.1167    6.0632   -0.9242   -0.4147
%     4.3058    5.1588    5.1466    5.1110    6.6811    3.9690    6.9090    5.3127    9.0561    6.6167    3.4116
%     7.9082    6.4769    7.2934    6.4178    5.3347    6.0745    7.8238    7.1962    6.1255    6.1808    4.2141
%     0.3993    2.0731    1.8289    2.5400    0.3058    1.5307    2.6737    1.9569    0.0034    2.1793   -0.5780
%     1.5183    1.4867    1.3963    1.5903    1.9142    1.5661    1.9727    2.0129    1.7221    1.7757    1.0986
%    -0.0991    0.1150   -1.1115    0.3519    1.5179   -0.6649    1.0695   -0.5581   -0.2070   -0.9709    0.6931];

ferikeNoScale200 = [  198.0755  197.2116  198.1600  197.2116  211.0693  197.5662  196.9025  202.6267  353.9279
    0.6503    0.1730    1.5187    0.1730   -0.1229    5.0164    0.8771    1.0165   -0.5386
    5.7534    6.1776    5.3462    6.1776    5.7932    6.0511    5.9253    5.3342    4.1183
    6.2000    5.5392    5.5254    5.5392    8.6764    5.5861    5.7226    5.9335    4.5138
    7.5431    5.7541    6.7479    5.7541    5.4606    5.6581    6.4203    5.9551    4.1037
    0.6648    0.9988    3.4932    0.9988    1.7476    0.3931    1.5579    2.0333   -0.4840
    5.2343    4.8669    5.7322    4.8669    5.2513    5.3185    4.9926    5.2592    3.3469
    4.8643    6.1236    4.8963    6.1236    5.6796    5.8015    5.5548    5.4513    4.2099
    0.9701    0.8060    0.5200    0.8060    2.9351    0.7145    1.1851   -1.4807   -0.5936
    1.3353    1.5166    1.5791    1.5166    1.5827    1.6938    1.7761    1.6410    1.0986
   -0.0162   -0.5427    0.3269   -0.5427    0.4467   -0.1364    0.1687   -0.8815    0.6931];

% ferikeNoScaleCentralized = [  124.7307  132.0423  127.4957  123.0596  129.4415  127.8795  129.5309  133.9958  123.0528  134.3121  140.6848
%     0.9154    0.8778   -0.0420    3.0416    0.2567    0.8779    1.0641   -0.0523   -0.2705   -3.1474   -0.4938
%     6.5921    5.6762    5.2440    7.2898    5.9144    5.0138    7.8011    5.9890    5.1785    7.2914    4.2097
%     5.2278    5.4092    5.2941    6.5143    5.8314    5.6994    5.5366    6.8321    5.8165    5.2352    4.5829
%     4.3783    7.4874    4.6593    6.3377    6.1940    4.9927    5.8131    5.0683    5.7835    2.1745    4.1260
%    -0.1242    0.9371    0.7819   -1.9665    0.1302   -2.1828   -0.6225   -0.3825    1.8148   -2.4719   -0.4147
%     5.0696    4.7216    4.9524    5.8570    4.9411    4.2111    4.8693    4.3952    6.7540   -0.2378    3.4116
%     6.3838    4.4833    5.6385    5.2525    5.6033    7.4610    5.2224    2.9917    5.2816    5.1484    4.2141
%     0.6238   -2.4891    0.3389    0.2656   -0.4258    1.3645    0.7368   -0.4052   -3.1474    2.3608   -0.5780
%     0.8172    0.6337    0.7167    0.9944    0.6916    1.0309    1.0986    0.7466    1.0451    0.7483    1.0986
%    -1.4109   -0.3143   -0.2773   -1.0723   -0.8282   -1.2893   -0.8065   -0.7879   -1.5640   -0.6423    0.6931];

% ferikeLast50 = [   56.6257   52.3989   50.4470   70.8683   50.7962   56.0267   59.3065   56.1461   94.8638
%     1.0478    0.4908    0.9050    0.7818    1.7133    1.7784    2.4760    0.1711   -0.4657
%     5.7910    7.1302    7.0687    6.2087    6.8230    6.0166    5.8400    5.9896    4.2573
%     5.7757    6.6915    6.6238    5.8243    6.3065    6.7789    6.7574    5.4214    4.5100
%     7.9401    6.4348    6.9273    4.0094    7.2845    8.3724    6.7799    7.0027    4.1021
%     0.9484    0.8964    1.3226   -0.2045    2.5182    0.2826    1.1285    1.5374   -0.5204
%     6.2655    7.6594    5.6530    3.7405    5.5348    5.3546    3.3010    8.3858    3.4080
%     7.4940    9.6124    9.0075    6.1840    8.3120    6.7014    6.2794    8.3877    4.2340
%     2.1239    0.3693    0.0431   -0.5557   -0.9239    0.4845    1.3011    1.6583   -0.6063
%     2.0489    1.7511    1.4438    1.8538    1.8163    1.7643    1.7433    1.5481    1.0986
%    -0.5031   -0.9487   -1.1620    0.4087   -0.4553    0.3009   -0.3637   -0.4562    0.6931];
foptsave = ferikeNoScale200(1, :);
lhoptsave = ferikeNoScale200(2:end, :)';
err = [];
% [~, ic] = sort(foptsave);
for i= 5%1:size(ferikeNoScale200, 2)
    
    
    loghyp = [lhoptsave(i, :)'; log(siga)];
    
    hyp = exp(loghyp);
    
    w = hyp(1:end-3); W = diag(w.^-2);
    sigf = hyp(end-2);
    sigp = hyp(end-1);
    siga = hyp(end);
    
    I = eye(size(samples, 1));
    Sigma = sigf^2 * exp(-.5 * maha([samples, context], [samples, context], W)) + I * ridge;
    
    activations = sum(Sigma/sigf^2, 2)/size(Sigma, 2);
    disp(['mean activation: ', num2str(mean(activations))])
    
    
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(zeros(size(samples, 1), 1), prefs, Sigma, sigp, fPrior, siga);
    
    fmapTarget = fmap;
    
    for j = 1:5
        ixTest = j:3:size(samples, 1);
        ixTrain = setdiff(1:size(samples, 1), ixTest);
        
        %     ixLoc = bsxfun(@eq, fPrior(:, 1), ixTrain);
        %
        %     I = eye(length(ixTrain));
        %     Sigma = sigf^2 * exp(-.5 * maha([samples(ixTrain, :), context(ixTrain, :)], [samples(ixTrain, :), context(ixTrain, :)], W)) + I * ridge;
        %
        %     activations = sum(Sigma/sigf^2, 2)/size(Sigma, 2);
        %     disp(['mean activation: ', num2str(mean(activations))])
        %
        %
        %     [fmap, ddS, GammaMap] = nr_plgp_wPrior(zeros(length(ixTrain), 1), prefs(ixLocPref, :), Sigma, sigp, fPrior(ixLocPrior, :), siga);
        
        % Sigma = I/(ddS - GammaMap);
        % Var = I/(Sigma + I/GammaMap);
        
        kall = sigf^2 *exp(-.5 * maha([samples(ixTrain, :) context(ixTrain, :)], [samples(ixTest, :) context(ixTest, :)], W));
        Kxx = sigf^2 *exp(-.5 * maha([samples(ixTest, :) context(ixTest, :)], [samples(ixTest, :) context(ixTest, :)], W)) ;
        %         SigmaStar = Kxx - kall' / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall;
        ypred = (kall' / Sigma(ixTrain, ixTrain)) * (fmap(ixTrain, :) );
        
        err(i,j) = mse(ypred - fmap(ixTest, :));
        
    end
    err2(i) = mse(fPrior(:, 2) - fmap(fPrior(:, 1), :));
end

% [~, ic] = sort(fmap);
combined = [[1:size(samples, 1)]', samples, context, fmap];
% combined = flipud(combined(ic, :));

combined = [combined, zeros(size(combined, 1), 5)];
for i = 1:size(forceJerk, 1)
    combined((combined(:, 1) == forceJerk(i, 1)), end-4:end) = forceJerk(i, 2:end);
end

ix1 = combined(:, 9) == 1;
ix1wm = and(ix1, combined(:, 12) ~=0);
ix2 = combined(:, 9) == 2;
ix2wm = and(ix2, combined(:, 12) ~=0);
ix3 = combined(:, 9) == 3;
ix3wm = and(ix3, combined(:, 12) ~=0);
ix4 = combined(:, 9) == 4;
ix4wm = and(ix4, combined(:, 12) ~=0);

ixStaticWm = (ix1wm + ix2wm)==1;
ixDynamicWm =  (ix3wm + ix4wm)==1;
ixAllWm =  (ix1wm +ix2wm +ix3wm + ix4wm)==1;


figure,
for i = 2:8
    subplot(2, 4, i)
    plot(combined(ix2, i), combined(ix2, 10), '*')
end

scaler = 3;
feri = [ combined(ix1, 2:8)' * exp(-(10-combined(ix1, 10)) * scaler)/sum(exp(-(10-combined(ix1, 10))* scaler)) ...
  combined(ix2, 2:8)' * exp(-(10-combined(ix2, 10))* scaler)/sum(exp(-(10-combined(ix2, 10))* scaler)) ...
   combined(ix3, 2:8)' * exp(-(10-combined(ix3, 10))* scaler)/sum(exp(-(10-combined(ix3, 10))* scaler)) ...
    combined(ix4, 2:8)' * exp(-(10-combined(ix4, 10))* scaler)/sum(exp(-(10-combined(ix4, 10))* scaler)) ...
    combined([ix1+ix2] ==1, 2:8)' * exp(-(10-combined([ix1+ix2] ==1, 10))* scaler)/sum(exp(-(10-combined([ix1+ix2] ==1, 10))* scaler)) ...
    combined([ix4+ix3] ==1, 2:8)' * exp(-(10-combined([ix4+ix3] ==1, 10))* scaler)/sum(exp(-(10-combined([ix4+ix3] ==1, 10))* scaler))]

% combinedOriginal = combined;
% combinedOriginal(:, 2:8) = bsxfun(@plus, bsxfun(@times, samples, sampleStd), sampleMean);
% combinedOriginal(:, 9)  = context .* contextStd + contextMean;
% combinedOriginal(:, 10) = combined(:, 10) + 5;
%
% relevances = zeros(8, 5);
%
%     loghyp = getFullGPModel(combined(ix1, 2:8), combined(ix1, 10), -200);
%     hyp = exp(loghyp); w= hyp(1:end-2);
%     relevances(1:7, 1) = w.^-2;
%
%     loghyp = getFullGPModel(combined(ix2, 2:8), combined(ix2, 10), -200);
%     hyp = exp(loghyp); w= hyp(1:end-2);
%     relevances(1:7, 2) = w.^-2 ;
%
%     loghyp = getFullGPModel(combined(ix3, 2:8), combined(ix3, 10), -200);
%     hyp = exp(loghyp); w= hyp(1:end-2);
%     relevances(1:7, 3) = w.^-2 ;
%
%     loghyp = getFullGPModel(combined(ix4, 2:8), combined(ix4, 10), -200);
%     hyp = exp(loghyp); w= hyp(1:end-2);
%     relevances(1:7, 4) = w.^-2 ;
%
%     loghyp = getFullGPModel(combined(:, 2:9), combined(:, 10), -200);
%     hyp = exp(loghyp); w= hyp(1:end-2);
%     relevances(:, 5) = w.^-2

figure, plot3(combined(ixDynamicWm, 14), combined(ixDynamicWm, 15), combined(ixDynamicWm, 10), '*'), grid on, xlabel('forcePeakTIme'), ylabel('handoverTime'), zlabel('reward'), title('dynamic handovers')
figure, plot(combined(ixDynamicWm, 15)-combined(ixDynamicWm, 14), combined(ixDynamicWm, 10), '*'), grid on, xlabel('firstPearkToHandover TIme'), ylabel('reward'), title('dynamic handovers')
figure, plot(combined(ixDynamicWm, 15), combined(ixDynamicWm, 10), '*'), grid on, xlabel('handover Time'), ylabel('reward'), title('dynamic handovers')

feri = [combined(ixDynamicWm, 15), combined(ixDynamicWm, 10)];
feri(:, 1) = feri(:, 1) + randn(size(feri, 1), 1)*.01;
[~, ic] = sort(feri(:, 1));
feri(ic, :)

xq = .4:.01:2;
jani = interp1(feri(:, 1), feri(:, 2), xq);

figure, plot3(combined(ixStaticWm, 14), combined(ixStaticWm, 15), combined(ixStaticWm, 10), '*'), grid on, xlabel('forcePeakTIme'), ylabel('handoverTime'), zlabel('reward'), title('static handovers')
figure, plot(combined(ixStaticWm, 15)-combined(ixStaticWm, 14), combined(ixStaticWm, 10), '*'), grid on, xlabel('firstPearkToHandover TIme'), ylabel('reward'), title('statichandovers')
figure, plot(combined(ixStaticWm, 15), combined(ixStaticWm, 10), '*'), grid on, xlabel('handover TIme'), ylabel('reward'), title('statichandovers')

figure, plot3(combined(ixAllWm, 14), combined(ixAllWm, 15), combined(ixAllWm, 10), '*'), grid on, xlabel('forcePeakTIme'), ylabel('handoverTime'), zlabel('reward'), title('allhandovers')
figure, plot(combined(ixAllWm, 15)-combined(ixAllWm, 14), combined(ixAllWm, 10), '*'), grid on, xlabel('firstPearkToHandover TIme'), ylabel('reward'), title('all handovers')
figure, plot(combined(ixAllWm, 15), combined(ixAllWm, 10), '*'), grid on, xlabel('handover TIme'), ylabel('reward'), title('all handovers')




[n, x] = hist(combined(ixAllWm, 15));
figure, plot(x, smooth(n/sum(n)), 'b')
hold on, 
ixAllWm = find(ixAllWm);
highRew = combined(ixAllWm, 10) >= mean(combined(ixAllWm, 10) );
[n, x] = hist(combined(ixAllWm(highRew), 15));
plot(x, smooth(n/sum(n)), 'r')

lowRew = combined(ixAllWm, 10) <= mean(combined(ixAllWm, 10) );
[n, x] = hist(combined(ixAllWm(lowRew), 15));
plot(x, smooth(n/sum(n)), 'k')








