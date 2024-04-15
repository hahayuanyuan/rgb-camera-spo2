param=[];
param.fps = 60;
param.L = param.fps*10;
param.bandwidth = [0.7 3];

filenam = [];
% filenam.video = '/Volumes/WD_BLACK/2303082040/video_cam.bin';
% filenam.video = '/Volumes/WD_BLACK/3th NICU/230302/2023hos3nicu13_230302/video_cam.bin';
filenam.video = '/Volumes/WD_BLACK/202303231230/video_cam.bin';
length = readBinFrameIds(filenam.video, 0);
% width = size(videoData,2);
% height = size(videoData,1);

trace = zeros(3,length);



% ROI selection

% 
%     Method 1: adjust ROI for each video
figure
rggb = readBinFrameIds(filenam.video, 600);
rgb = mydebayer(rggb);
im = rgb;
imagesc(uint8(im(:,:,:)*8));
roi = drawpolygon('Color', [0 0.4470 0.7410]);
mask = createMask(roi);



for index=1:length
    rggb = readBinFrameIds(filenam.video, index);
    tmpImage = mydebayer(rggb);
%     tmpImage = videoData(:,:,:,index);

    r = tmpImage(:,:,1);
    g = tmpImage(:,:,2);
    b = tmpImage(:,:,3);

    
    m_r = mean(r(mask == 1));
    m_g = mean(g(mask == 1));
    m_b = mean(b(mask == 1));
    trace(:,index) = [m_r; m_g; m_b];
end
%%
trace = trace(:,1:14000);

figure
subplot(3,1,1),plot(trace(1,:),'Color',[0.8 0.2 0.2],'LineWidth',1.5);
set(gca,'FontSize',15)
subplot(3,1,2),plot(trace(2,:),'Color',[0.2 0.8 0.2],'LineWidth',1.5);
set(gca,'FontSize',15)
subplot(3,1,3),plot(trace(3,:),'Color',[0.2 0.2 0.8],'LineWidth',1.5);
set(gca,'FontSize',15)

log4trace = log(trace);




r1 = zeros(1,size(log4trace,2)-param.L+1);
r2 = zeros(1,size(log4trace,2)-param.L+1);
r3 = zeros(1,size(log4trace,2)-param.L+1);
r4 = zeros(1,size(log4trace,2)-param.L+1);
rr1 = zeros(2,size(log4trace,2)-param.L+1);
rr2 = zeros(2,size(log4trace,2)-param.L+1);
rr3 = zeros(2,size(log4trace,2)-param.L+1);
rr4 = zeros(2,size(log4trace,2)-param.L+1);
rr5 = zeros(1,size(log4trace,2)-param.L+1);
for ni = 1:size(log4trace,2)-param.L+1
    partialC = log4trace(:,ni:ni+param.L-1);
        
    [n, d] = butter(2, param.bandwidth/(param.fps/2), 'bandpass');
    partialCf = filtfilt(n,d,partialC')';
    partialCf = partialCf .* repmat(hann(size(partialCf,2))',[3,1]);
    
    [rr,r] = acAcquired(partialCf);
    rr1(:,ni) = rr(:,1);
    rr2(:,ni) = rr(:,2);
    rr3(:,ni) = rr(:,3);
    rr4(:,ni) = rr(:,4);
    r1(ni) = r(1);
    r2(ni) = r(2);
    r3(ni) = r(3);
    r4(ni) = r(4);

    if ni <= 1000
        rr5(ni) = rr(1,4);
    else
%         amplitude4r = abs(fft(partialCf(1,:)));
%         amplitude4g = abs(fft(partialCf(2,:)));
%         ratio4test = amplitude4r./amplitude4g;
%         error4test = abs(ratio4test-repmat(median(rr5(ni-60:ni-1)),[1,param.L]));
%         [~,index4test] = min(error4test);
%         rr5(ni) = ratio4test(index4test);

    end

end

% figure
% subplot(3,1,1),plot(log4trace(1,:));
% subplot(3,1,2),plot(log4trace(2,:));
% subplot(3,1,3),plot(log4trace(3,:));

cf4log4trace = filtfilt(n,d,log4trace')';

% figure
% subplot(3,1,1),plot(cf4log4trace(1,:));
% subplot(3,1,2),plot(cf4log4trace(2,:));
% subplot(3,1,3),plot(cf4log4trace(3,:));


figure
subplot(4,1,1)
plot(param.L/2:param.L/2+size(rr1,2)-1,rr1(1,:),'Color','k','LineWidth',1.5)
% plot(201:2900,rr1(1,201:2900),'Color','k','LineWidth',1.5)
title('Median Peak - Median Valley')
set(gca,'FontSize',15)

subplot(4,1,2)
plot(param.L/2:param.L/2+size(rr2,2)-1,rr2(1,:),'Color','k','LineWidth',1.5)
% plot(201:2900,rr2(1,201:2900),'Color','k','LineWidth',1.5)
title('Square Power')
set(gca,'FontSize',15)

subplot(4,1,3)
plot(param.L/2:param.L/2+size(rr3,2)-1,rr3(1,:),'Color','k','LineWidth',1.5)
% plot(201:2900,rr3(1,201:2900),'Color','k','LineWidth',1.5)
title('Shannon Power')
set(gca,'FontSize',15)

subplot(4,1,4)
plot(param.L/2:param.L/2+size(rr4,2)-1,rr4(1,:),'Color','k','LineWidth',1.5)
% plot(201:2900,rr4(1,201:2900),'Color','k','LineWidth',1.5)
title('Hilbert')
set(gca,'FontSize',15)

figure
plot(rr5)
% figure
% subplot(4,1,1)
% % plot(log4trace(1,:))
% % hold on
% plot(param.L/2:param.L/2+size(r1,2)-1,r1)
% 
% subplot(4,1,2)
% % plot(log4trace(1,:))
% % hold on
% plot(param.L/2:param.L/2+size(r2,2)-1,r2)
% 
% subplot(4,1,3)
% % plot(log4trace(1,:))
% % hold on
% plot(param.L/2:param.L/2+size(r3,2)-1,r3)
% 
% subplot(4,1,4)
% plot(param.L/2:param.L/2+size(r4,2)-1,r4)

function [ratio, r, g, b] = acAcquired(partialCf)


    gap=20;
    [maxR, ~] = findpeaks(partialCf(1,:),'MinPeakDistance',gap);
    [maxG, ~] = findpeaks(partialCf(2,:),'MinPeakDistance',gap);
    [maxB, ~] = findpeaks(partialCf(3,:),'MinPeakDistance',gap);


    [minR, ~] = findpeaks(-partialCf(1,:),'MinPeakDistance',gap);
    [minG, ~] = findpeaks(-partialCf(2,:),'MinPeakDistance',gap);
    [minB, ~] = findpeaks(-partialCf(3,:),'MinPeakDistance',gap);
    minR = -minR;
    minG = -minG;
    minB = -minB;
    

    varianceR1 = median(maxR)-median(minR);
    varianceG1 = median(maxG)-median(minG);
    varianceB1 = median(maxB)-median(minB);
    

    ratio1 = [varianceR1/varianceG1;varianceR1/varianceB1];

    % normalized to [-1,1]
    maxr = max(abs(partialCf(1,:)));
    maxg = max(abs(partialCf(2,:)));
    maxb = max(abs(partialCf(3,:)));
%     maxrgb = max([max(abs(partialCf(1,:))),max(abs(partialCf(2,:))),max(abs(partialCf(3,:)))]);
    partialCf(1,:) = partialCf(1,:)./maxr;
    partialCf(2,:) = partialCf(2,:)./maxg;
    partialCf(3,:) = partialCf(3,:)./maxb;

    tmpr = partialCf(1,:).*partialCf(1,:);
    tmpg = partialCf(2,:).*partialCf(2,:);
    tmpb = partialCf(3,:).*partialCf(3,:);

    [n1, d1] = butter(1, 0.1/(60/2), 'low');
    acr2 = filtfilt(n1,d1,tmpr')';
    acg2 = filtfilt(n1,d1,tmpg')';
    acb2 = filtfilt(n1,d1,tmpb')';

    varianceR2 = mean(acr2)*maxr;
    varianceG2 = mean(acg2)*maxg;
    varianceB2 = mean(acb2)*maxb;

    ratio2 = [varianceR2/varianceG2;varianceR2/varianceB2];

    tmpr = tmpr+1e-9;
    tmpg = tmpg+1e-9;
    tmpb = tmpb+1e-9;

    tmpr = -tmpr.*log(tmpr);
    tmpg = -tmpg.*log(tmpg);
    tmpb = -tmpb.*log(tmpb);

    acr3 = filtfilt(n1,d1,tmpr')';
    acg3 = filtfilt(n1,d1,tmpg')';
    acb3 = filtfilt(n1,d1,tmpb')';

    varianceR3 = mean(acr3)*maxr;
    varianceG3 = mean(acg3)*maxg;
    varianceB3 = mean(acb3)*maxb;
    
    ratio3 = [varianceR3/varianceG3;varianceR3/varianceB3];

    zr = hilbert(partialCf(1,:));
    zg = hilbert(partialCf(2,:));
    zb = hilbert(partialCf(3,:));

    varianceR4 = median(abs(zr))*maxr;
    varianceG4 = median(abs(zg))*maxg;
    varianceB4 = median(abs(zb))*maxb;

    ratio4 = [varianceR4/varianceG4;varianceR4/varianceB4];


    ratio = [ratio1,ratio2,ratio3,ratio4];
    r = [varianceR1,varianceR2,varianceR3,varianceR4];

    
end
