function [i] = SaveMovie(u, x, xmin, xmax, maximum, t, filename, nplots, i)
% Save movie of u to file with filename handle

% Form frames and save
figure(i)
clear M;
video = VideoWriter(filename);
video.FrameRate = 15;
for j = 1:nplots
    plot(x,real(u(j,:)),x,imag(u(j,:)),x,abs(u(j,:)));
    text(pi/2,3*maximum/4,strcat('t = ',num2str(t(j))))
    
    % Window depending on choice of initial data
    axis([xmin xmax -0.3*maximum 1.2*maximum]),
    M(j) = getframe(gcf);  
end
open(video);
writeVideo(video,M);
close(video);

% Increment i
i = i+1;


end

