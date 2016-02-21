% bidirectionally align 3D movie
function movieL = bidi_align(movie)
temp = zeros(size(movie,3),1);

reduction = size(movie,2)/8;
cross_template = zeros(size(movie,2)*2-1-4*reduction,1)';
cross_template2 = zeros(size(movie,2)*2-1-4*reduction,1)';
for j = 1:size(movie,3)
    j/size(movie,3)
    for o = 1:size(movie,1)/2
        dd1 = movie((o-1)*2+1,:,j); dd1 = dd1(reduction+1:end-reduction);
        dd2 = movie((o-1)*2+2,:,j); dd2 = dd2(reduction+1:end-reduction);
        cross_template = cross_template + xcorr(dd1,dd2,'biased');
        cross_template2 = cross_template2 + xcorr(dd1,dd1,'biased');
    end
    [~,a] = max(cross_template);
    [~,b] = max(cross_template2);
    temp(j) = a-b; 
end

movieL = zeros(size(movie));

for j = 1:size(movie,3)
    for o = 1:size(movie,1)/2
        movieL((o-1)*2+2,:,j) = circshift(movie((o-1)*2+2,:,j),[0 0]);
        movieL((o-1)*2+1,:,j) = circshift(movie((o-1)*2+1,:,j),[0 -mean(temp)]);
    end
end
disp(strcat('Mean change:',12,num2str(mean(temp)),12,'Std change:',12,num2str(std(temp))));

end