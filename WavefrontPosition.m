function [ wavefront, wavefrontindex, i] = ... 
    WavefrontPosition( u, N, t, x, xmin, xmax, relativecutoff, method, ... 
    wavenumber, i, amp0, width, toplot)
% Calculate approximate location of wavefront position in moving front DNLS
% solution with N grid points using either spatial or spectral method with
% prescribed cutoff value

% Define midpoint
midpoint = (xmax - xmin) / 2;

if strcmp(method,'spatial')
    
    % Initialize wavefront vector
    wavefront = zeros(1,length(t));
    wavefront(:) = midpoint + width;
    wavefrontindex = zeros(1,length(t));
    for j = 1:length(t)
        l = N/2;
        while l <= N
            if abs(u(j,l))/min(abs(u(j,:))) >= relativecutoff
                l = l+1;
            else
                wavefront(j) = x(l);
                wavefrontindex(j) = l;
                l = N+1;
            end
        end
        wavefrontindex(j) = max(wavefrontindex(j),1);
    end
end

if strcmp(method,'spectral')
    
    wavefront = zeros(1,length(t));
    wavefront(:) = midpoint + width;
    wavefrontindex = zeros(1,length(t));
    phase = angle(u);

    for l = 1: length(t)
        wavenumber(l,:) = deriv((phase(l,:)'),N,1)';
    end

    for j = 1:length(t)
        l = N/2;
        while l <= N
            if abs(wavenumber(j,l))/mean(abs(wavenumber(j,:))) <= relativecutoff
                l = l+1;
            else
                wavefront(j) = x(l);
                wavefrontindex(j) = l;
                l = N+1;
            end
        end
        wavefrontindex(j) = max(wavefrontindex(j),1);
    end
end

% Plot wavefront position if requested
if strcmp(toplot,'yes')
    figure(i), plot(t,wavefront), 
    title(strcat('Wavefront Position: k = ', int2str(log2(N)), ... 
        ' , amplitude = ', num2str(amp0), ' , width =  ', num2str(width), ... 
        ' , method: ', method)),
    xlabel('t'), ylabel('Wavefront position')

    % Increment i
    i = i+1;
end

end

