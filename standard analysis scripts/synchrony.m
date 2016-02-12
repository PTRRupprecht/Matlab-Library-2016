
% global mode & synchrony

sFF = FF_0X(:,trash);

for kf = 1:4
fix = find( (sMI == kf)' & ~isnan(sum(SM,1)));


Data = FF_0X(:,trash);
Data = Data(:,fix);
% hh = find(~isnan(sum(Data)));
% Data = Data(:,hh);
% Data_shuffled = FF_0X(:,trash);
w_size = 500;
w_shift = 50;
nb_bins = (floor(size(Data,1)/w_shift)-w_size/w_shift);
clear chi2 chi
for j = 1:nb_bins
    j
    range = (1:w_size) + (j-1)*w_shift;
    rangeX(j) = round(mean(range));
    chi2(1) = NaN;
    for N = 2:size(Data,2)
        temp = 0;
        for k = 1:floor(numel(trash)/N)
            indizes = randi(size(Data,2),N,1);
            L = Data(range,indizes);
%             for u = 1:numel(indizes)
%                 shift = randi(numel(range),1);
%                 L(:,u) = circshift(L(:,u),[shift 0]);
%             end
            Vt = nanmean(L,2);
            sigma_v_2 = nanmean((Vt.^2)) - nanmean(Vt)^2;
            sigma_v_2i = nanmean(L.^2,1) - nanmean(L,1).^2;
            temp = temp + sigma_v_2/nanmean(sigma_v_2i)/ floor(numel(trash)/N);
        end
        chi2(N) = temp;
    end
    chi(j,:) = sqrt(chi2);
end

cmap = jet(nb_bins);
figure(88), hold on;
for j = 1:nb_bins
    warning off;
%     plot(log(1:numel(trash)),smooth(log(chi(j,:)),7),'Color',cmap(j,:));
    plot((1:size(Data,2)),smooth(chi(j,:),1),'Color',cmap(j,:));
    ft = fittype( 'a/sqrt(x) + b', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    % opts.Display = 'Off';
    % opts.Lower = [-Inf -Inf];
    % opts.StartPoint = [0.775177503995417 0.945009650034365];
    % opts.Upper = [Inf Inf];
    try
    % Fit model to data.
    [xData, yData] = prepareCurveData( 1:size(Data,2), chi(j,:) );

    [fitresult, gof] = fit(xData, yData, ft, opts );
    res1(j) = fitresult.a;
    res2(j) = fitresult.b;
    catch 
        res1(j) = NaN;
        res2(j) = NaN;
    end
end
hold off;

% figure, plot(rangeX,res2)

figure(kf*112), plot(rangeX,res2)

end

