%%======================================================
%                         EXO 2
% ======================================================
clear all;
close all;
clc;

%% Question 1

% Image =imread('omni.png');
% NB = rgb2gray(Image);
% 
% figure
% subplot 211
% imshow(Image)
% subplot 212
% imshow(NB)
% [x0, y0] = ginput;
% 
% [X2,Y2] = meshgrid(1:360,1:200);
% X1 = Y2.*cosd(X2)+x0;
% Y1 = Y2.*sind(X2)+y0;
% 
% J = interp2(double(NB),X1,Y1);
% 
% R = interp2(double(Image(:,:,1)),X1,Y1);
% G = interp2(double(Image(:,:,2)),X1,Y1);
% B = interp2(double(Image(:,:,3)),X1,Y1);
% img=cat(3,R,G,B);
% 
% figure
% subplot 211
% imshow(uint8(J))
% subplot 212
% imshow(uint8(img))

%% =====================================================
% Question 2
% 1)

% I=double(imread('lena256.png'));
% s = size(I);
% H = [1 0 10;0 1 20;0 0 1];
% Img1 = transformimage(I, H, s);
% 
% % 2)
% 
% x0=100 ; 
% y0=100 ;
% r = 0.4 ; % facteur d??chelle
% HT = [1 0 x0; 0 1 y0; 0 0 1];
% HS = [r 0 0 ; 0 r 0 ; 0 0 1];
% H = HT * HS * inv(HT);
% 
% Img2 = transformimage(I, H, s);
% 
% % 3)
% 
% x0=100 ; y0=100 ;
% a = pi/180* 20 ; % 20 degr?s
% HT = [1 0 x0; 0 1 y0; 0 0 1];
% HR = [cos(a) sin(a) 0 ; -sin(a) cos(a) 0 ; 0 0 1];
% H = HT * HR * inv(HT);
% 
% Img3 = transformimage(I, H, s);

% figure
% imshow(uint8(I))
% 
% figure
% imshow(uint8(Img1))
% 
% figure
% imshow(uint8(Img2))
% 
% figure
% imshow(uint8(Img3))

%% ====================================================
% Question 3
% a)
% [Ih, Iw]=size(I);
% pts1=[1 1; Iw 1; 1 Ih; Iw Ih];
% figure(1); clf
% imagesc(I);
% title('Image initiale');
% hold on
% plot(pts1(:,1),pts1(:,2),'o-r','linewidth',2);
% 
% % b)
% Jh=256; Jw=256;
% figure(2); clf
% imagesc(zeros(Jh,Jw));
% [x,y]=ginput(4);
% pts2=[x(:),y(:)];
% hold on
% plot(pts2(:,1),pts2(:,2),'o-b','linewidth',2);
% 
% % c)
% T=cp2tform(pts1, pts2,'projective');
% H=T.tdata.T';
% 
% J=transformimage(I, H, [Jh Jw]);
% 
% figure
% imshow(uint8(J))

%% ====================================================
% Question Bonus
% I = double(imread('qr-code-pub.jpg'));
% Im = double(imread('lena256.png'));
% 
% figure(1)
% imagesc(uint8(I));
% title('Image initiale');
% hold on
% 
% % b)
% Jh=256; Jw=256;
% x = [67.8 226.6 56 236];
% y = [60.2 29.8 219.5 210.8]
% pts1=[x(:),y(:)];
% pts2=[1 1; Jw 1; 1 Jh; Jw Jh];
% hold on
% plot(pts1(:,1),pts1(:,2),'o-b','linewidth',2);
% 
% % c)
% T=cp2tform(pts1, pts2,'projective');
% H=T.tdata.T';
% 
% J=transformimage(I, H, [Jh Jw]);
% 
% figure
% imshow(uint8(J))

%% ==================================================

I = double(imread('qr-code-pub.jpg'));
Im = double(imread('lena256.png'));

% b)
Jh=256; Jw=256;
x = [67.8 226.6 56 236];
y = [60.2 29.8 219.5 210.8]
pts1=[x(:),y(:)];
pts2=[1 1; Jw 1; 1 Jh; Jw Jh];

% c)
T=cp2tform(pts2, pts1,'projective');
H=T.tdata.T';

J=transformimage(Im, H, [Jh Jw]);
J = cat(3, J, J, J);

%mask = find(J == 0);
mask = zeros(Jh, Jw);

for i = 1:Jh
    for j = 1:Jw
        if J(i,j,1) == 0 & J(i,j,2) == 0 & J(i,j,3) == 0
            mask(i,j) = 1;
        end
        if mask(i,j) == 1
            J(i,j,:) = I(i,j,:);
        end
    end
end

figure
imshow(uint8(I));
hold on
plot(pts1(:,1),pts1(:,2),'o-b','linewidth',2);

figure
imshow(uint8(I))
hold on
imshow(uint8(J))
