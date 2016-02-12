

% load mp3s

FileList = dir('*.mp3');

[y,Fs] = audioread(FileList(6).name);



% downsampling

yy = y;
for k = 1:7
    yy(k:7:end,1) = y(1:7:end-k+1,1)+y(1:7:end-k+1,1);
    yy(k:7:end,2) = y(1:7:end-k+1,2);
end

yy(:,1) = smooth(yy(:,1),10);
yy(:,2) = smooth(yy(:,2),10);

LUT = round(linspace(500,1,numel(1:10:(size(y,1)-10000))));
counter = 1;
for i = 500:10:(size(y,1)-1000)
    y(i:(i+LUT(counter)*2),1) = smooth(y(i:(i+LUT(counter)*2),1),LUT(counter));
    y(i:(i+LUT(counter)*2),2) = smooth(y(i:(i+LUT(counter)*2),2),LUT(counter));
    counter = counter + 1;
end
    

y = (y + circshift(y,[1 0]) + circshift(y,[2 0]) + circshift(y,[3 0]) + circshift(y,[4 0]))/5;

yy = y(1:5:end,:);

audiowrite('rewritten3.wav',yy,Fs,'BitsPerSample',24);


% spectrogram

N = 1024;
n = 0:N-1;

w0 = 2*pi/5;
x = sin(w0*n)+10*sin(2*w0*n.^6);
figure(2), plot(x)
s = spectrogram(x);

figure(1); spectrogram(y(100000:440000,1),44000)

