function track_peak(int_cqt, x1, y1, x2,y2,w)

figure; hold on;

m = (y2-y1)/2;
m = round(y1+m);
accv = [(m-w):(m+w)]*0; 
 
for i=x1:x2;
    spec = int_cqt(i,(m-w):(m+w));
    %plot(spec);
    accv = accv+spec;
end
plot(accv);
figure;

