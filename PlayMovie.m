function [ i ] = PlayMovie( u, x, xmin, xmax, maximum, t, frames, loops, i )
% Play movie of solution to PDE that loops repeatedly

figure(i)
clear M;
% Fill frames
for j = 1:loops*frames
    if mod(j,frames) ~= 0
    plot(x,real(u(mod(j,frames),:)),x,imag(u(mod(j,frames),:)), ...
        x,abs(u(mod(j,frames),:))), 
    text(pi/2,3*maximum/4,strcat('t = ',num2str(t(mod(j,frames)))))
    axis([xmin xmax -0.3*maximum 1.2*maximum])
    M(j) = getframe(gcf);
    else
    plot(x,real(u(mod(j,frames)+frames,:)), ...
        x, imag(u(mod(j,frames)+frames,:)), x, ...
        abs(u(mod(j,frames)+frames,:))),
    text(pi/2,3*maximum/4,strcat('t = ',num2str(t(mod(j,frames)+frames))))
    axis([xmin xmax -0.3*maximum 1.2*maximum])
	M(j) = getframe(gcf);  
    end
end

% Play movie
movie(M);

% Increment i
i = i+1;

end

