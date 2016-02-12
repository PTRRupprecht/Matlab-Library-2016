% bidirectionally align 3D movie
function movieL = bidi_align(movie)
temp = zeros(size(movie,3),1);

cross_template = zeros(size(movie,2)*2-1,1)';
cross_template2 = zeros(size(movie,2)*2-1,1)';
for j = 1:size(movie,3)
    j/size(movie,3)
    for o = 1:size(movie,1)/2
        cross_template = cross_template + xcorr(movie((o-1)*2+1,:,j),movie((o-1)*2+2,:,j),'biased');
        cross_template2 = cross_template2 + xcorr(movie((o-1)*2+1,:,j),movie((o-1)*2+1,:,j),'biased');
    end
    [~,a] = max(cross_template);
    [~,b] = max(cross_template2);
    temp(j) = a-b; 
end

movieL = zeros(size(movie));

for j = 1:size(movie,3)
    for o = 1:size(movie,1)/2
        movieL((o-1)*2+2,:,j) = circshift(movie((o-1)*2+2,:,j),[0 0]);
        movieL((o-1)*2+1,:,j) = circshift(movie((o-1)*2+1,:,j),[0 b-a]);
    end
end
disp(strcat('Mean change:',12,num2str(mean(temp)),12,'Std change:',12,num2str(std(temp))));

end