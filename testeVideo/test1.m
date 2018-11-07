%# figure
%figure, set(gcf, 'Color','white')
%Z = peaks; surf(Z);  axis tight
%set(gca, 'nextplot','replacechildren', 'Visible','off');

%# create AVI object
nFrames = 20;
vidObj = VideoWriter('myPeaks3.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 10;
open(vidObj);

t = 0:0.01:2*pi;
%# create movie
for k=1:nFrames
   h= figure;
   plot(t,sin(2*pi*k*t))
   title(['Frequency: ' num2str(k)])
   xlabel('t')
   ylabel('sin')
   writeVideo(vidObj, getframe(h));
   close 
end
close(gcf)

%# save as AVI file, and open it using system video player
close(vidObj);