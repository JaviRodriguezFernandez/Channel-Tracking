% Problem 2 Gaussian Sums Filter (GSF)

clc, clear all, close all;

rng(1);
mx_0 = zeros(4,1);
Pbar_0 = diag([1,1e-6*ones(1,3)]);

time = 1:200;
NMC = 100;
N = 100;

x0 = mx_0 + chol(Pbar_0)*randn(4,1);

x_MC = zeros(4,length(time),NMC);
x_hat_MC = zeros(4,length(time),NMC);
P_hat_MC = zeros(4,4,length(time),NMC);
RMSE = zeros(4,NMC);

for nmc = 1:NMC
    Percentage_MC = nmc/NMC*100
    x0 = mx_0 + chol(Pbar_0)*randn(4,1);

    [x,x_hat,P_hat] = RBPF(time,x0,N);

    x_MC(:,:,nmc) = x;
    x_hat_MC(:,:,nmc) = x_hat;
    P_hat_MC(:,:,:,nmc) = P_hat;
    
    RMSE(:,nmc) = sqrt(mean(abs(x_hat-x).^2,2));

    % Evaluate Gaussian kernel density to find PDF evolution with time
    
%     figure(1), subplot(221), hold on,
%     plot(time,x(1,:),'*-',time,x_hat(1,:),'o-')
% %     plot(time,squeeze(sqrt(P_hat(1,1,:))),'d-')
%     hold off;
%     legend('True state x(1)','Estimate x(1)','Standard deviation')
%     
%     subplot(222), hold on,
%     plot(time,x(2,:),'*-',time,x_hat(2,:),'o-')
% %     plot(time,squeeze(sqrt(P_hat(2,2,:))),'d-')
%     hold off;
%     legend('True state x(2)','Estimate x(2)','Standard deviation')
%     
%     subplot(223), hold on,
%     plot(time,x(3,:),'*-',time,x_hat(3,:),'o-')
% %     plot(time,squeeze(sqrt(P_hat(3,3,:))),'d-')
%     hold off;
%     legend('True state x(3)','Estimate x(3)','Standard deviation')
%     
%     subplot(224), hold on,
%     plot(time,x(4,:),'*-',time,x_hat(4,:),'o-')
% %     plot(time,squeeze(sqrt(P_hat(4,4,:))),'d-')
%     hold off;
%     legend('True state x(4)','Estimate x(4)','Standard deviation')
end

RMSE_av = mean(RMSE,2)


x_MC_1 = squeeze(x_MC(:,end,:));
mean_x = squeeze(mean(x_MC_1,2));
mean_x_1 = mean_x(1:2);
mean_x_2 = mean_x(3:4);

sample_cov = zeros(4,4);

for nmc = 1:NMC
    sample_cov = sample_cov + (x_MC_1(:,nmc) - mean_x)*(x_MC_1(:,nmc) - mean_x).';
end
sample_cov_x = 1/NMC*sample_cov;
sample_cov_x = diag(sample_cov_x);

% PDF for real state
    
Sampling_density = 100;

Support_density_orig = mean_x_1 + [linspace(-3*sqrt(sample_cov_x(1)),3*sqrt(sample_cov_x(1)),Sampling_density); linspace(-3*sqrt(sample_cov_x(2)),3*sqrt(sample_cov_x(2)),Sampling_density)];
Support_density = [repelem(Support_density_orig(1,:),Sampling_density); repmat(Support_density_orig(2,:),1,Sampling_density)].';
Sample_PDF = ksdensity(x_MC_1(1:2,:).',Support_density);
Grid_Sample_PDF_1 = reshape(Sample_PDF,Sampling_density,Sampling_density);
figure, subplot(211), surf5 = surf(Support_density_orig(1,:),Support_density_orig(2,:),Grid_Sample_PDF_1,'LineStyle','none')
title('Estimated joint PDF of x_1 and x_2 at the final time');
Support_density_orig = mean_x_2 + [linspace(-3*sqrt(sample_cov_x(3)),3*sqrt(sample_cov_x(3)),Sampling_density); linspace(-3*sqrt(sample_cov_x(4)),3*sqrt(sample_cov_x(4)),Sampling_density)];
Support_density = [repelem(Support_density_orig(1,:),Sampling_density); repmat(Support_density_orig(2,:),1,Sampling_density)].';
Sample_PDF = ksdensity(x_MC_1(3:4,:).',Support_density);
Grid_Sample_PDF_2 = reshape(Sample_PDF,Sampling_density,Sampling_density);
subplot(212), surf6 = surf(Support_density_orig(1,:),Support_density_orig(2,:),Grid_Sample_PDF_2,'LineStyle','none')
title('Estimated joint PDF of x_3 and x_4 at the final time');

% PDF for estimate state


x_hat_MC_1 = squeeze(x_hat_MC(:,end,:));
mean_x = squeeze(mean(x_hat_MC_1,2));
mean_x_1 = mean_x(1:2);
mean_x_2 = mean_x(3:4);

sample_cov = zeros(4,4);

for nmc = 1:NMC
    sample_cov = sample_cov + (x_hat_MC_1(:,nmc) - mean_x)*(x_hat_MC_1(:,nmc) - mean_x).';
end
sample_cov_x = 1/NMC*sample_cov;
sample_cov_x = diag(sample_cov_x);
    
Sampling_density = 100;

Support_density_orig = mean_x_1 + [linspace(-3*sqrt(sample_cov_x(1)),3*sqrt(sample_cov_x(1)),Sampling_density); linspace(-3*sqrt(sample_cov_x(2)),3*sqrt(sample_cov_x(2)),Sampling_density)];
Support_density = [repelem(Support_density_orig(1,:),Sampling_density); repmat(Support_density_orig(2,:),1,Sampling_density)].';
Sample_PDF = ksdensity(x_MC_1(1:2,:).',Support_density);
Grid_Sample_PDF_1 = reshape(Sample_PDF,Sampling_density,Sampling_density);
figure, subplot(211), surf5 = surf(Support_density_orig(1,:),Support_density_orig(2,:),Grid_Sample_PDF_1,'LineStyle','none')
title('Estimated joint PDF of \hat{x}_1 and \hat{x}_2 at the final time');

Support_density_orig = mean_x_2 + [linspace(-3*sqrt(sample_cov_x(3)),3*sqrt(sample_cov_x(3)),Sampling_density); linspace(-3*sqrt(sample_cov_x(4)),3*sqrt(sample_cov_x(4)),Sampling_density)];
Support_density = [repelem(Support_density_orig(1,:),Sampling_density); repmat(Support_density_orig(2,:),1,Sampling_density)].';
Sample_PDF = ksdensity(x_MC_1(3:4,:).',Support_density);
Grid_Sample_PDF_2 = reshape(Sample_PDF,Sampling_density,Sampling_density);
subplot(212), surf6 = surf(Support_density_orig(1,:),Support_density_orig(2,:),Grid_Sample_PDF_2,'LineStyle','none')
title('Estimated joint PDF of \hat{x}_3 and \hat{x}_4 at the final time');




Error = x_MC - x_hat_MC;

sample_mean = 1/NMC*sum(Error,3);
sample_cov = zeros(4,4,length(time));
for t = 1:length(time)
    for nmc = 1:NMC
        sample_cov(:,:,t) = sample_cov(:,:,t) + 1/NMC*(Error(:,t,nmc) - sample_mean(:,t))*(Error(:,t,nmc) - sample_mean(:,t)).';
    end
end

figure(1), subplot(221), hold on,
plot1(1:NMC) = plot(time,squeeze(Error(1,:,:)),'b-')
plot2(1:NMC) = plot(time,3*squeeze(sqrt(P_hat_MC(1,1,:,:))),'r-')
plot3(1:NMC) = plot(time,-3*squeeze(sqrt(P_hat_MC(1,1,:,:))),'r-')
plot4(1) = plot(time,squeeze(sample_mean(1,:)),'k*-')
plot5(1) = plot(time,3*squeeze(sqrt(sample_cov(1,1,:))),'m*-')
plot6(1) = plot(time,-3*squeeze(sqrt(sample_cov(1,1,:))),'m*-')
legend([plot1(1),plot2(1),plot4(1),plot5(1)],{'Estimation error','\pm 3 x Predicted Standard Deviation','Sample mean','\pm 3 x Sample standard deviation'});
xlabel('t')
ylabel('Statistics of estimation error')
title('Component # 1')
axis tight
set(gca,...
  'Box'         , 'on',...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...  
    'LineWidth',2,...
    'FontSize'    , 14);
grid on;

hold off;

subplot(222), hold on,
plot1(1:NMC) = plot(time,squeeze(Error(2,:,:)),'b-')
plot2(1:NMC) = plot(time,3*squeeze(sqrt(P_hat_MC(2,2,:,:))),'r-')
plot3(1:NMC) = plot(time,-3*squeeze(sqrt(P_hat_MC(2,2,:,:))),'r-')
plot4(1) = plot(time,squeeze(sample_mean(2,:)),'k*-')
plot5(1) = plot(time,3*squeeze(sqrt(sample_cov(2,2,:))),'m*-')
plot6(1) = plot(time,-3*squeeze(sqrt(sample_cov(2,2,:))),'m*-')
legend([plot1(1),plot2(1),plot4(1),plot5(1)],{'Estimation error','\pm 3 x Predicted Standard Deviation','Sample mean','\pm 3 x Sample standard deviation'});
xlabel('t')
ylabel('Statistics of estimation error')
title('Component # 2')
axis tight
set(gca,...
  'Box'         , 'on',...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...  
    'LineWidth',2,...
    'FontSize'    , 14);
grid on;

hold off;

subplot(223), hold on,
plot1(1:NMC) = plot(time,squeeze(Error(3,:,:)),'b-')
plot2(1:NMC) = plot(time,3*squeeze(sqrt(P_hat_MC(3,3,:,:))),'r-')
plot3(1:NMC) = plot(time,-3*squeeze(sqrt(P_hat_MC(3,3,:,:))),'r-')
plot4(1) = plot(time,squeeze(sample_mean(3,:)),'k*-')
plot5(1) = plot(time,3*squeeze(sqrt(sample_cov(3,3,:))),'m*-')
plot6(1) = plot(time,-3*squeeze(sqrt(sample_cov(3,3,:))),'m*-')
legend([plot1(1),plot2(1),plot4(1),plot5(1)],{'Estimation error','\pm 3 x Predicted Standard Deviation','Sample mean','\pm 3 x Sample standard deviation'});
xlabel('t')
ylabel('Statistics of estimation error')
title('Component # 3')
axis tight
set(gca,...
  'Box'         , 'on',...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...  
    'LineWidth',2,...
    'FontSize'    , 14);
grid on;

hold off;

subplot(224), hold on,
plot1(1:NMC) = plot(time,squeeze(Error(4,:,:)),'b-')
plot2(1:NMC) = plot(time,3*squeeze(sqrt(P_hat_MC(4,4,:,:))),'r-')
plot3(1:NMC) = plot(time,-3*squeeze(sqrt(P_hat_MC(4,4,:,:))),'r-')
plot4(1) = plot(time,squeeze(sample_mean(4,:)),'k*-')
plot5(1) = plot(time,3*squeeze(sqrt(sample_cov(4,4,:))),'m*-')
plot6(1) = plot(time,-3*squeeze(sqrt(sample_cov(4,4,:))),'m*-')
legend([plot1(1),plot2(1),plot4(1),plot5(1)],{'Estimation error','\pm 3 x Predicted Standard Deviation','Sample mean','\pm 3 x Sample standard deviation'});
xlabel('t')
ylabel('Statistics of estimation error')
title('Component # 4')
axis tight
set(gca,...
  'Box'         , 'on',...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...  
    'LineWidth',2,...
    'FontSize'    , 14);
grid on;

hold off;