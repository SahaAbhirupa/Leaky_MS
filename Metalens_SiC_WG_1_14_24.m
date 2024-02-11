% Metalense Phase Profile
%clear all

% 633 nm lambda
% 100 micron focal length
l = 633e-9;
f = 100e-6;
[x,y] = meshgrid(-40:0.315:40,-40:0.624:40);
% y= meshgrid(-40:0.415:40,-40:0.415:40);
[m, n] = size(x) ;

n_eff=2.025;
% % phi=importdata(phase_x.csv);
% Degrees
p = mod((360/l)*(f-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2))-((n_eff*360*x*10^-6/l)),360);
% p = mod((360/l)*(-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)),360);
for i=1:size(x,1)
    for j=1:size(x,2)
        if(p(i,j)>180)
            p(i,j) = p(i,j)-360;
            % p(i,j)= p(i,j)- phi(i,j);
        end
    end
end

%phi = -1i*sind(p);
%Radians
%p = mod((-2*pi/l)*(f-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)),2*pi);
%for i=1:length(x)
%    for j=1:length(x)
%        if(p(i,j)>pi)
%            p(i,j) = p(i,j)-(2*pi);
%        end
%    end
%end
%p = -1i*p;

% Without modulus 
%p = -1i*(-2*pi/l)*(f-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2));

% Calculate final phase portion of E-field
phi = exp(real(p)).*(cos(imag(p))+i*sin(imag(p)));
figure(1)
subplot(1,2,1)
imagesc(p)
colormap('default')
colorbar
title('Metalens Phase Profile')
xlabel('# of meta-atom')
ylabel('# of meta-atom')
view(2)
pbaspect([1 1 1])
% Fourier Transform
% Angular Spectrum
P = (fftshift(fft2(phi))).^2;
P = abs(P/max(P(:)));
%figure(2)
subplot(1,2,2)

imagesc(P)
title('Intensity Profile')
colormap('default')
colorbar

pbaspect([1 1 1])
figure(3)
%subplot(3,1,3)

plot(P(65,:))
title('Intensity Profile Cross Section')
xlabel('# of meta-atom')
ylabel('Normalized Intensity')
xlim([0 160])
pbaspect([1 1 1])
% Inverse Fourier Transform
%IP = (ifft2(P));
%figure(3)
%imagesc(real(IP))
%colormap('default')
%colorbar
        
% Remove side lobes
% option 1: airy phase and grating filter
%p = mod(((x+y).^3+(y-x).^3)*(0.5),360);
%[x,y] = meshgrid(-20:0.5:20,-20:0.1:20);
%[m, n] = size(x);
p = mod((360/l)*(f-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2))-((n_eff*360*x*10^-6/l)),360);

for i=1:size(x,1)
    for j=1:size(x,2)
        if(p(i,j)>180)
            p(i,j) = p(i,j)-360;
        end
    end
end
writematrix(p,'phase_sic_633nm_try.txt');
writematrix(p,'phase_sic_633nm_try.csv');
figure(4)

imagesc(p)
colormap('jet')
colorbar
title('Metalens Phase Profile')
xlabel('# of meta-atom')
ylabel('# of meta-atom')
view(2)
phi = exp(real(p)).*(cos(imag(p))+i*sin(imag(p)));
P = (fftshift(fft2(phi))).^2;
figure(5)
imagesc(abs(P/max(P(:))))
title('Intensity Profile')
colormap('default')
colorbar

%Adding Linear phase grating
lim = 30; %where the filter starts
scale = 360; % filter scale
%[x,y] = meshgrid(-20:0.5:20,-20:0.5:20);
%[m, n] = size(x);
% p = mod((360/l)*(-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)-(n_eff*360*x/l)),360);

% p = mod((360/l)*(-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)+(n_eff*360*x*10^-6/l)),360);
p = mod((360/l)*(f-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)),360);
for i=1:length(x)
    for j=1:size(x,2)
        if abs(x(i,j))>=lim && abs(y(i,j))<=lim
            p(i,j) = mod((y(i,j)*scale),360);
        elseif abs(y(i,j))>=lim && abs(x(i,j))<=lim
            p(i,j) = mod((x(i,j)*scale),360);  
        elseif abs(y(i,j))>=lim && abs(x(i,j))>=lim
            if sign(x(i,j))==sign(y(i,j))
                p(i,j) = mod((((x(i,j))-y(i,j))*scale),360);
            else
                p(i,j) = mod((((x(i,j))+y(i,j))*scale),360);
            end
        end
        if(p(i,j)>180)
            p(i,j) = p(i,j)-360;
        end
    end
end
figure(6)
subplot(1,2,1)
imagesc(p)
colormap('default')
colorbar
title('Metalens Phase Profile')
xlabel('# of meta-atom')
ylabel('# of meta-atom')
view(2)
pbaspect([1 1 1])
phi = exp(real(p)).*(cos(imag(p))+i*sin(imag(p)));
P = (fftshift(fft2(phi)))^2;
P = abs(P/max(P(:)));
%figure(7)
subplot(1,2,2)
imagesc(P)
title('Intensity Profile')
colormap('default')
colorbar     
pbaspect([1 1 1])
figure(8)
plot(P(81,:))
title('Intensity Profile Cross Section')
xlabel('# of meta-atom')
ylabel('Normalized Intensity')     
xlim([0 160])
% Test out more filters
lim = 35; %where the filter starts
scale = 360/2; % filter scale
% p = mod((360/l)*(-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)-(n_eff*360*x*10^-6/l)),360);
p = mod((360/l)*(-sqrt((x*10^-6).^2+(y*10^-6).^2+f^2)),360);
for i=1:length(x)
    for j=1:length(x)
        dist = sqrt(x(i,j).^2+y(i,j).^2);
        if dist >= lim
            p(i,j) = 72*(mod(atand((y(i,j)./x(i,j))),5)-2.5);
        end
        if(p(i,j)>180)
            p(i,j) = p(i,j)-360;
        end
    end
end
figure(9)
imagesc(p)
writematrix(p,'phase_sic_633nm_corrected.txt');
colormap('default')
colorbar
title('Metalens Phase Profile')
xlabel('# of meta-atom')
ylabel('# of meta-atom')
view(2)
phi = exp(real(p)).*(cos(imag(p))+i*sin(imag(p)));
P = (fftshift(fft2(phi)))^2;
figure(10)
imagesc(abs(P/max(P(:))))
title('Intensity Profile')
colormap('default')
colorbar             
        
        