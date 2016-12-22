function [ rx, phix, i ] = MoreSnapshots( r, phi, x, t, num_snapshots, N, i )
% Plot snapshots of other functions of modulus and phase at different
% timesteps

% Instantiate rx and phix variables
rx = zeros(size(r));
phix = zeros(size(phi));

% Calculate spatial derivatives of modulus and phase to fill rx, phix
for j = 1:length(t)
   phix(j,:) = deriv(phi(j,:)',N,1)'; 
end
for j = 1:length(t)
   rx(j,:) = deriv(r(j,:)',N,1)';
end

% Group velocity
gpvelocity = r.^2.*phix;

% r^2 r_x
r2rx = r.^2.*rx;

% Compute integral of r^3 phi_x^2 and compare to r^2r_x wavefront jump
othertermderiv = r.^3.*phix.^2;
otherterm = zeros(size(othertermderiv));
meshsize = x(2)-x(1);

% Numerically calculate integral from 0 to x of othertermderiv
for k = 1:length(t)
    currentothertermtrunc = zeros(N/2,1);
    for l = 1:N/2
        currenttestpoint = othertermderiv(k,N/2+l);
        if l == 1
            currentothertermtrunc(l) = 0;
        else
            currentothertermtrunc(l) = currentothertermtrunc(l-1)+ ... 
                currenttestpoint*meshsize;
        end
        
 
    end
    otherterm(k,:) = [wrev(currentothertermtrunc); ... 
    currentothertermtrunc];   
%     currentothertermderivtrunc = othertermderiv(k,N/2+1:N);
%     xtrunc = x(N/2+1:N);
%     [xnew,currentothertermtrunc] = ... 
%         ode23(@IntRHS,xtrunc,currentothertermderivtrunc(1));
%     otherterm(k,:) = ... 
%         [wrev(currentothertermtrunc);currentothertermtrunc];
end



% Introduce number of snapshots
increment = (length(t)-1)/num_snapshots;

% Come up with subplot arrangement
verplots = floor(sqrt(num_snapshots));
horplots = floor(num_snapshots/verplots);
while horplots*verplots ~= num_snapshots
    verplots = verplots + 1;
    horplots = floor(num_snapshots/verplots);
end


for j = 1:horplots
    for k = 1:verplots
        % Find current time
        currenttime = (((verplots*(j-1))+k-1) * increment)+1;
        
        % Plot r^2 phi_x snapshot
        figure(i)
        subplot(horplots,verplots,verplots*(j-1)+k)
        plot(x,abs(gpvelocity(currenttime,:)))
        xlim([0 2*pi])
        title(['r^2 \phi_x: t = ', num2str(t(currenttime))])
        
        % Plot r^2 r_x snapshot
        figure(i+1)
        subplot(horplots,verplots,verplots*(j-1)+k)
        plot(x,abs(r2rx(currenttime,:)))
        xlim([0 2*pi])
        title(['r^2 r_x: t = ', num2str(t(currenttime))])
        
        % Plot r^3 phi_x^2 snapshot
        figure(i+2)
        subplot(horplots,verplots,verplots*(j-1)+k)
        plot(x,abs(othertermderiv(currenttime,:)))
        xlim([0 2*pi])
        title(['r^3 \phi_x^2: t = ', num2str(t(currenttime))])
        
        % Plot integral of r^3 phi_x^2 snapshot
        figure(i+3)
        subplot(horplots,verplots,verplots*(j-1)+k)
        plot(x,abs(otherterm(currenttime,:)))
        xlim([0 2*pi])
        title(['\int r^3 \phi_x^2: t = ', num2str(t(currenttime))])
    end
end

% Increment i
i = i+4;


end

