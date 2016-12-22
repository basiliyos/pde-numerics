function [r, phi, i] = ... 
    Snapshots( u, x, xmin, xmax, maximum, t, num_snapshots, cutoff, i )
% Plot modulus and phase of solution at different timesteps

% Create and fill r, phi vectors
indexvec = 1:floor((length(t))/num_snapshots):length(t);
r = abs(u);
phi = zeros(size(u));
for j = 1:length(t)
    for l = 1:length(x)
        if r(j,l) > cutoff
            phi(j,l) = angle(u(j,l));
        end
    end
end

% Fill and plot r and phi at various times
figure(i)
for j = 1:num_snapshots
    
    % Plot modulus snapshots
    subplot(2,num_snapshots,j)
    plot(x,r(indexvec(j),:))
    axis([xmin, xmax,-0.3*maximum,1.2*maximum])
    title(['Modulus: t = ', num2str(t(indexvec(j)))])
     
    % Plot phase snapshots
    subplot(2,num_snapshots,j+5)
    plot(x,phi(indexvec(j),:))
    axis([xmin, xmax, -pi, pi])
    title(['Phase: t = ', num2str(t(indexvec(j)))])
end

% Increment i
i = i+1;

end

