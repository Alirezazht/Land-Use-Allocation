clc;clear all; close all;
ax = imread('asli.tif');

maskoni   = find(ax==1);
tejari    = find(ax==2);
amozeshi  = find(ax==3);
khadamati = find(ax==4);
sanati    = find(ax==5);
park      = find(ax==6);
keshavarzi= find(ax==7);
bayer     = find(ax==8);
khali     = find(ax==0);
nonkhali  = find(ax~=0);

karbari = find(ax~=0);

karbari2 =zeros(length(karbari),2);
karbari2(:,1)=karbari;
karbari2(:,2)=ax(karbari);
tedad={zeros(8,8)};
tedad=repmat(tedad,[256 249]);

for i=1:256
    for j=1:249
        for ring=1:8
            for x=i-ring:i+ring
                if x<=0 || x>256
                    continue
                end
                for y=j-ring:j+ring
                    if y<=0 || y>249
                        continue
                    end
                    if x==i-ring || x==i+ring || y==j-ring || y==j+ring
                        if ax(x,y)==1
                            tedad{i,j}(1,ring)= tedad{i,j}(1,ring)+1;
                        elseif ax(x,y)==2
                            tedad{i,j}(2,ring)= tedad{i,j}(2,ring)+1;
                        elseif ax(x,y)==3
                            tedad{i,j}(3,ring)= tedad{i,j}(3,ring)+1;
                        elseif ax(x,y)==4
                            tedad{i,j}(4,ring)= tedad{i,j}(4,ring)+1;
                        elseif ax(x,y)==5
                            tedad{i,j}(5,ring)= tedad{i,j}(5,ring)+1;
                        elseif ax(x,y)==6
                            tedad{i,j}(6,ring)= tedad{i,j}(6,ring)+1;
                        elseif ax(x,y)==7
                            tedad{i,j}(7,ring)= tedad{i,j}(7,ring)+1;
                        elseif ax(x,y)==8
                            tedad{i,j}(8,ring)= tedad{i,j}(8,ring)+1;     
                        end
                    end
                end
            end
        end
    end
end


%MASKONI
M_karbari = karbari2(karbari2(:,2)==1,1);
F_M={zeros(6,8)};
F_M=repmat(F_M,[1,length(M_karbari)]);
F_M_M=zeros(8,8);
r_pix=8:8:64;
for i=1:length(M_karbari)
    for ring=1:8
        F_M{i}(1,ring)=(tedad{M_karbari(i)}(1,ring)/r_pix(ring))/(length(maskoni)/32595);
        F_M{i}(2,ring)=(tedad{M_karbari(i)}(2,ring)/r_pix(ring))/(length(tejari)/32595);
        F_M{i}(3,ring)=(tedad{M_karbari(i)}(3,ring)/r_pix(ring))/(length(amozeshi)/32595);
        F_M{i}(4,ring)=(tedad{M_karbari(i)}(4,ring)/r_pix(ring))/(length(khadamati)/32595);
        F_M{i}(5,ring)=(tedad{M_karbari(i)}(5,ring)/r_pix(ring))/(length(sanati)/32595);
        F_M{i}(6,ring)=(tedad{M_karbari(i)}(6,ring)/r_pix(ring))/(length(park)/32595);
        F_M{i}(7,ring)=(tedad{M_karbari(i)}(7,ring)/r_pix(ring))/(length(keshavarzi)/32595);
        F_M{i}(8,ring)=(tedad{M_karbari(i)}(8,ring)/r_pix(ring))/(length(bayer)/32595);
    end
    F_M_M=F_M_M+F_M{i};
end
F_M_M=F_M_M/length(M_karbari);
F_M_M1=log10(F_M_M);
figure
plot(1:8,F_M_M1(1,:))
grid on
title('maskoni bar maskoni')
figure
plot(1:8,F_M_M1(2,:))
grid on
title('tejari bar maskoni')
figure
plot(1:8,F_M_M1(3,:))
grid on
title('amozeshi bar maskoni')
figure
plot(1:8,F_M_M1(4,:))
grid on
title('khadamati bar maskoni')
figure
plot(1:8,F_M_M1(5,:))
grid on
title('sanaati bar maskoni')
figure
plot(1:8,F_M_M1(6,:))
grid on
title('park bar maskoni')
figure
plot(1:8,F_M_M1(7,:))
grid on
title('keshavarzi bar maskoni')
figure
plot(1:8,F_M_M1(8,:))
grid on
title('bayer bar maskoni')

%TEJARI
T_karbari=karbari2(karbari2(:,2)==2,1);
F_T={zeros(6,8)};
F_T=repmat(F_T,[1,length(T_karbari)]);
F_T_M=zeros(8,8);
r_pix=8:8:64;
for i=1:length(T_karbari)
    for ring=1:8
        F_T{i}(1,ring)=(tedad{T_karbari(i)}(1,ring)/r_pix(ring))/(length(maskoni)/32595);
        F_T{i}(2,ring)=(tedad{T_karbari(i)}(2,ring)/r_pix(ring))/(length(tejari)/32595);
        F_T{i}(3,ring)=(tedad{T_karbari(i)}(3,ring)/r_pix(ring))/(length(amozeshi)/32595);
        F_T{i}(4,ring)=(tedad{T_karbari(i)}(4,ring)/r_pix(ring))/(length(khadamati)/32595);
        F_T{i}(5,ring)=(tedad{T_karbari(i)}(5,ring)/r_pix(ring))/(length(sanati)/32595);
        F_T{i}(6,ring)=(tedad{T_karbari(i)}(6,ring)/r_pix(ring))/(length(park)/32595);
        F_T{i}(7,ring)=(tedad{T_karbari(i)}(7,ring)/r_pix(ring))/(length(keshavarzi)/32595);
        F_T{i}(8,ring)=(tedad{T_karbari(i)}(8,ring)/r_pix(ring))/(length(bayer)/32595);
    end
    F_T_M=F_T_M+F_T{i};
end
F_T_M=F_T_M/length(T_karbari);
F_T_M1=log10(F_T_M);
figure
plot(1:8,F_T_M1(1,:))
grid on
title('maskoni bar tejari')
figure
plot(1:8,F_T_M1(2,:))
grid on
title('tejari bar tejari')
figure
plot(1:8,F_T_M1(3,:))
grid on
title('amozeshi bar tejari')
figure
plot(1:8,F_T_M1(4,:))
grid on
title('khadamati bar tejari')
figure
plot(1:8,F_T_M1(5,:))
grid on
title('sanaati bar tejari')
figure
plot(1:8,F_T_M1(6,:))
grid on
title('park bar tejari')
figure
plot(1:8,F_T_M1(7,:))
grid on
title('keshavarzi bar tejari')
figure
plot(1:8,F_T_M1(8,:))
grid on
title('bayer bar tejari')

%AMOZESHI
A_karbari =karbari2(karbari2(:,2)==3,1);
F_A={zeros(6,8)};
F_A=repmat(F_A,[1,length(A_karbari)]);
F_A_M=zeros(8,8);
r_pix=8:8:64;
for i=1:length(A_karbari)
    for ring=1:8
        F_A{i}(1,ring)=(tedad{A_karbari(i)}(1,ring)/r_pix(ring))/(length(maskoni)/32595);
        F_A{i}(2,ring)=(tedad{A_karbari(i)}(2,ring)/r_pix(ring))/(length(tejari)/32595);
        F_A{i}(3,ring)=(tedad{A_karbari(i)}(3,ring)/r_pix(ring))/(length(amozeshi)/32595);
        F_A{i}(4,ring)=(tedad{A_karbari(i)}(4,ring)/r_pix(ring))/(length(khadamati)/32595);
        F_A{i}(5,ring)=(tedad{A_karbari(i)}(5,ring)/r_pix(ring))/(length(sanati)/32595);
        F_A{i}(6,ring)=(tedad{A_karbari(i)}(6,ring)/r_pix(ring))/(length(park)/32595);
        F_A{i}(7,ring)=(tedad{A_karbari(i)}(7,ring)/r_pix(ring))/(length(keshavarzi)/32595);
        F_A{i}(8,ring)=(tedad{A_karbari(i)}(8,ring)/r_pix(ring))/(length(bayer)/32595);
    end
    F_A_M=F_A_M+F_A{i};
end
F_A_M=F_A_M/length(A_karbari);
F_A_M1=log10(F_A_M);
figure
plot(1:8,F_A_M1(1,:))
grid on
title('maskoni bar amozeshi')
figure
plot(1:8,F_A_M1(2,:))
grid on
title('tejari bar amozeshi')
figure
plot(1:8,F_A_M1(3,:))
grid on
title('amozeshi bar amozeshi')
figure
plot(1:8,F_A_M1(4,:))
grid on
title('khadamati bar amozeshi')
figure
plot(1:8,F_A_M1(5,:))
grid on
title('sanaati bar amozeshi')
figure
plot(1:8,F_A_M1(6,:))
grid on
title('park bar amozeshi')
figure
plot(1:8,F_A_M1(7,:))
grid on
title('keshavarzi bar amozeshi')
figure
plot(1:8,F_A_M1(8,:))
grid on
title('bayer bar amozeshi')



%KHADAMATI
K_karbari=karbari2(karbari2(:,2)==4,1);
F_K={zeros(6,8)};
F_K=repmat(F_K,[1,length(K_karbari)]);
F_K_M=zeros(8,8);
r_pix=8:8:64;
for i=1:length(K_karbari)
    for ring=1:8
        F_K{i}(1,ring)=(tedad{K_karbari(i)}(1,ring)/r_pix(ring))/(length(maskoni)/32595);
        F_K{i}(2,ring)=(tedad{K_karbari(i)}(2,ring)/r_pix(ring))/(length(tejari)/32595);
        F_K{i}(3,ring)=(tedad{K_karbari(i)}(3,ring)/r_pix(ring))/(length(amozeshi)/32595);
        F_K{i}(4,ring)=(tedad{K_karbari(i)}(4,ring)/r_pix(ring))/(length(khadamati)/32595);
        F_K{i}(5,ring)=(tedad{K_karbari(i)}(5,ring)/r_pix(ring))/(length(sanati)/32595);
        F_K{i}(6,ring)=(tedad{K_karbari(i)}(6,ring)/r_pix(ring))/(length(park)/32595);
        F_K{i}(7,ring)=(tedad{K_karbari(i)}(7,ring)/r_pix(ring))/(length(keshavarzi)/32595);
        F_K{i}(8,ring)=(tedad{K_karbari(i)}(8,ring)/r_pix(ring))/(length(bayer)/32595);
    end
    F_K_M=F_K_M+F_K{i};
end
F_K_M=F_K_M/length(K_karbari);
F_K_M1=log10(F_K_M);
figure
plot(1:8,F_K_M1(1,:))
grid on
title('maskoni bar khadamati')
figure
plot(1:8,F_K_M1(2,:))
grid on
title('tejari bar khadamati')
figure
plot(1:8,F_K_M1(3,:))
grid on
title('amozeshi bar khadamati')
figure
plot(1:8,F_K_M1(4,:))
grid on
title('khadamati bar khadamati')
figure
plot(1:8,F_K_M1(5,:))
grid on
title('sanaati bar khadamati')
figure
plot(1:8,F_K_M1(6,:))
grid on
title('park bar khadamati')
figure
plot(1:8,F_K_M1(7,:))
grid on
title('keshavarzi bar khadamati')
figure
plot(1:8,F_K_M1(8,:))
grid on
title('bayer bar khadamati')


%sanaati
S_karbari=karbari2(karbari2(:,2)==5,1);
F_S={zeros(6,8)};
F_S=repmat(F_S,[1,length(S_karbari)]);
F_S_M=zeros(8,8);
r_pix=8:8:64;
for i=1:length(S_karbari)
    for ring=1:8
        F_S{i}(1,ring)=(tedad{S_karbari(i)}(1,ring)/r_pix(ring))/(length(maskoni)/32595);
        F_S{i}(2,ring)=(tedad{S_karbari(i)}(2,ring)/r_pix(ring))/(length(tejari)/32595);
        F_S{i}(3,ring)=(tedad{S_karbari(i)}(3,ring)/r_pix(ring))/(length(amozeshi)/32595);
        F_S{i}(4,ring)=(tedad{S_karbari(i)}(4,ring)/r_pix(ring))/(length(khadamati)/32595);
        F_S{i}(5,ring)=(tedad{S_karbari(i)}(5,ring)/r_pix(ring))/(length(sanati)/32595);
        F_S{i}(6,ring)=(tedad{S_karbari(i)}(6,ring)/r_pix(ring))/(length(park)/32595);
        F_S{i}(7,ring)=(tedad{S_karbari(i)}(7,ring)/r_pix(ring))/(length(keshavarzi)/32595);
        F_S{i}(8,ring)=(tedad{S_karbari(i)}(8,ring)/r_pix(ring))/(length(bayer)/32595);
    end
    F_S_M=F_S_M+F_S{i};
end
F_S_M=F_S_M/length(S_karbari);
F_S_M1=log10(F_S_M);
figure
plot(1:8,F_S_M1(1,:))
grid on
title('maskoni bar sanaati')
figure
plot(1:8,F_S_M1(2,:))
grid on
title('tejari bar sanaati')
figure
plot(1:8,F_S_M1(3,:))
grid on
title('amozeshi bar sanaati')
figure
plot(1:8,F_S_M1(4,:))
grid on
title('khadamati bar sanaati')
figure
plot(1:8,F_S_M1(5,:))
grid on
title('sanaati bar sanaati')
figure
plot(1:8,F_S_M1(6,:))
grid on
title('park bar sanaati')
figure
plot(1:8,F_S_M1(7,:))
grid on
title('keshavarzi bar sanaati')
figure
plot(1:8,F_S_M1(8,:))
grid on
title('bayer bar sanaati')

%park
P_karbari=karbari2(karbari2(:,2)==6,1);
F_P={zeros(6,8)};
F_P=repmat(F_P,[1,length(P_karbari)]);
F_P_M=zeros(8,8);
r_pix=8:8:64;
for i=1:length(P_karbari)
    for ring=1:8
        F_P{i}(1,ring)=(tedad{P_karbari(i)}(1,ring)/r_pix(ring))/(length(maskoni)/32595);
        F_P{i}(2,ring)=(tedad{P_karbari(i)}(2,ring)/r_pix(ring))/(length(tejari)/32595);
        F_P{i}(3,ring)=(tedad{P_karbari(i)}(3,ring)/r_pix(ring))/(length(amozeshi)/32595);
        F_P{i}(4,ring)=(tedad{P_karbari(i)}(4,ring)/r_pix(ring))/(length(khadamati)/32595);
        F_P{i}(5,ring)=(tedad{P_karbari(i)}(5,ring)/r_pix(ring))/(length(sanati)/32595);
        F_P{i}(6,ring)=(tedad{P_karbari(i)}(6,ring)/r_pix(ring))/(length(park)/32595);
        F_P{i}(7,ring)=(tedad{P_karbari(i)}(7,ring)/r_pix(ring))/(length(keshavarzi)/32595);
        F_P{i}(8,ring)=(tedad{P_karbari(i)}(8,ring)/r_pix(ring))/(length(bayer)/32595);
    end
    F_P_M=F_P_M+F_P{i};
end
F_P_M=F_P_M/length(P_karbari);
F_P_M1=log10(F_P_M);

figure
plot(1:8,F_P_M1(1,:))
grid on
title('maskoni bar park')
figure
plot(1:8,F_P_M1(2,:))
grid on
title('tejari bar park')
figure
plot(1:8,F_P_M1(3,:))
grid on
title('amozeshi bar park')
figure
plot(1:8,F_P_M1(4,:))
grid on
title('khadamati bar park')
figure
plot(1:8,F_P_M1(5,:))
grid on
title('sanaati bar park')
figure
plot(1:8,F_P_M1(6,:))
grid on
title('park bar park')
figure
plot(1:8,F_S_M1(7,:))
grid on
title('keshavarzi bar park')
figure
plot(1:8,F_P_M1(8,:))
grid on
title('bayer bar park')

%-----------------------------------------------------
%-----------------------------------------------------
%-----------------------------------------------------
F={F_M_M1,
   F_T_M1,
   F_A_M1,
   F_K_M1,
   F_S_M1,
   F_P_M1};

for k=1:6
    for x=1:6
        for y=1:8
            if F{k}(x,y)==inf || F{k}(x,y)==-inf
                F{k}(x,y)=0;
            end
        end
    end
end

karbari3={zeros(256,249)};
karbari3=repmat(karbari3,[1,6]);

for L=1:6
    for i=1:256
        for j=1:249
            if ax(i,j)~=0

                for ring=1:8
                    for xx=i-ring:i+ring
                        if xx<=0 || xx>256
                            continue
                        end
                        for yy=j-ring:j+ring
                            if yy<=0 || yy>249
                                continue
                            end
                            if xx==i-ring || xx==i+ring || yy==j-ring || yy==j+ring
                                if ax(xx,yy)==1
                                    karbari3{L}(i,j)= karbari3{L}(i,j)+F{L}(1,ring);
                                elseif ax(xx,yy)==2
                                    karbari3{L}(i,j)= karbari3{L}(i,j)+F{L}(2,ring);
                                elseif ax(xx,yy)==3
                                    karbari3{L}(i,j)= karbari3{L}(i,j)+F{L}(3,ring);
                                elseif ax(xx,yy)==4
                                    karbari3{L}(i,j)= karbari3{L}(i,j)+F{L}(4,ring);
                                elseif ax(xx,yy)==5
                                    karbari3{L}(i,j)= karbari3{L}(i,j)+F{L}(5,ring);
                                elseif ax(xx,yy)==6
                                    karbari3{L}(i,j)= karbari3{L}(i,j)+F{L}(6,ring);
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


xxx =zeros(2,6);
karbari_normal={zeros(256,249)};
karbari_normal=repmat(karbari_normal,[1,6]);
for L=1:6
    Z=karbari3{L}(karbari);
    xxx(1,L)=(max(max(karbari3{L}))-min(min(Z)))/0.8;
    xxx(2,L)=(xxx(1,L)/10)-min(min(karbari3{L}));
    for i=karbari
        karbari_normal{L}(i)=(karbari3{L}(i)+xxx(2,L))/xxx(1,L);
        karbari_normal{L}(i)=min(max(karbari_normal{L}(i),0),1);
    end
end
% -------------------------------------------------------------------------
figure
imagesc(karbari_normal{1})
title('maskoni')
figure
imagesc(karbari_normal{2})
title('tejari')
figure
imagesc(karbari_normal{3})
title('amozeshi')
figure
imagesc(karbari_normal{4})
title('khadamati')
figure
imagesc(karbari_normal{5})
title('sanaati')
figure
imagesc(karbari_normal{6})
title('park')

imwrite(karbari_normal{1},'maskoni.tiff')
imwrite(karbari_normal{2},'tejari.tiff')
imwrite(karbari_normal{3},'amozeshi.tiff')
imwrite(karbari_normal{4},'khadamati.tiff')
imwrite(karbari_normal{5},'sanaati.tiff')
imwrite(karbari_normal{6},'park.tiff')