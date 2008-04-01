function[]=readfile(nameoffile1,nameoffile2,nameoffile3)

%clear all;
%reads info stored in files by dopper program
fid = fopen(nameoffile1, 'r');

%reading the variables from the file, in order
sizeofsearch=fscanf(fid, '%lg', [1 1]);

for i=1:sizeofsearch
RAsourcevector(i)=fscanf(fid, '%lg', [1 1]);
Decsourcevector(i)=fscanf(fid, '%lg', [1 1]);
powermatrix1(i)=fscanf(fid, '%lg', [1 1]);
end
RAsourcevector=RAsourcevector*360/2/pi;
Decsourcevector=Decsourcevector*360/2/pi;

frequency1= fscanf(fid, '%lg', [1 1])
%     amp = fscanf(fid, '%lg', [1 1])
meanofnoise= fscanf(fid, '%lg', [1 1])
SNR= fscanf(fid, '%lg', [1 1])
stdofnoise= fscanf(fid, '%lg', [1 1])
maximum= fscanf(fid, '%lg', [1 1])
rightascension_of_source= fscanf(fid, '%lg', [1 1])
declination_of_source= fscanf(fid, '%lg', [1 1])
RAsource=fscanf(fid, '%lg', [1 1])
Decsource=fscanf(fid, '%lg', [1 1])
TRI = delaunay(RAsourcevector,Decsourcevector);
nface=size(TRI);
for i=1:nface(1) 
color(i)=(1/3*(powermatrix1(TRI(i,1))+powermatrix1(TRI(i,2))+powermatrix1(TRI(i,3))));
color(i)=abs(color(i));
end
f1=figure(1)
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on','ParallelLabel','on','LabelFormat','none','MeridianLabel','on','MLabelParallel','equator');
patchm(Decsourcevector,RAsourcevector,'Faces',TRI,'FaceColor','flat','LineStyle','none','FaceVertexCData',color')
xlabel('Right Ascension');
ylabel('Declination');
title({'SNR OVER THE SKY',' freq=',frequency1,'Hz,  ', '   mean power of noise =',meanofnoise,'   SNR =', SNR, '   stdv of noise = ', stdofnoise,'maxpower =', maximum,'RA,DEC(deg) =', rightascension_of_source, declination_of_source  });
colorbar('location','southoutside');

plotm([declination_of_source],[rightascension_of_source],'x');
plotm([Decsource*360/2/pi],[RAsource*360/2/pi],'+');
fclose(fid);
saveas(f1,'f1.fig')



%clear all;
%reads info stored in files by dopper program
fid = fopen(nameoffile2, 'r');

%reading the variables from the file, in order
sizeofsearch=fscanf(fid, '%lg', [1 1]);

for i=1:sizeofsearch
RAsourcevector(i)=fscanf(fid, '%lg', [1 1]);
Decsourcevector(i)=fscanf(fid, '%lg', [1 1]);
powermatrix2(i)=fscanf(fid, '%lg', [1 1]);
end
RAsourcevector=RAsourcevector*360/2/pi;
Decsourcevector=Decsourcevector*360/2/pi;
frequency2= fscanf(fid, '%lg', [1 1])
%amp = fscanf(fid, '%lg', [1 1])
meanofnoise= fscanf(fid, '%lg', [1 1])
SNR= fscanf(fid, '%lg', [1 1])
stdofnoise= fscanf(fid, '%lg', [1 1])
maximum= fscanf(fid, '%lg', [1 1])
rightascension_of_source= fscanf(fid, '%lg', [1 1])
declination_of_source= fscanf(fid, '%lg', [1 1])
RAsource=fscanf(fid, '%lg', [1 1])
Decsource=fscanf(fid, '%lg', [1 1])
TRI = delaunay(RAsourcevector,Decsourcevector);
nface=size(TRI);
for i=1:nface(1) 
color(i)=(1/3*(powermatrix2(TRI(i,1))+powermatrix2(TRI(i,2))+powermatrix2(TRI(i,3))));
color(i)=abs(color(i));
end
f2=figure(2)
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on','ParallelLabel','on','LabelFormat','none','MeridianLabel','on','MLabelParallel','equator');patchm(Decsourcevector,RAsourcevector,'Faces',TRI,'FaceColor','flat','LineStyle','none','FaceVertexCData',color')
xlabel('Right Ascension');
ylabel('Declination');
title({'SNR OVER THE SKY',' freq=',frequency2, '   mean power of noise =',meanofnoise,'   SNR =', SNR, '   stdv of noise = ', stdofnoise,'maxpower =', maximum,'RA,DEC(deg) =', rightascension_of_source, declination_of_source  });
colorbar('location','southoutside');
plotm([declination_of_source],[rightascension_of_source],'x');
plotm([Decsource*360/2/pi],[RAsource*360/2/pi],'+');
fclose(fid);
saveas(f2,'f2.fig')




%clear all;
%reads info stored in files by dopper program
fid = fopen(nameoffile3, 'r');

%reading the variables from the file, in order
sizeofsearch=fscanf(fid, '%lg', [1 1]);

for i=1:sizeofsearch
RAsourcevector(i)=fscanf(fid, '%lg', [1 1]);
Decsourcevector(i)=fscanf(fid, '%lg', [1 1]);
powermatrix3(i)=fscanf(fid, '%lg', [1 1]);
end
RAsourcevector=RAsourcevector*360/2/pi;
Decsourcevector=Decsourcevector*360/2/pi;
frequency3= fscanf(fid, '%lg', [1 1])
%     amp = fscanf(fid, '%lg', [1 1])
meanofnoise= fscanf(fid, '%lg', [1 1])
SNR= fscanf(fid, '%lg', [1 1])
stdofnoise= fscanf(fid, '%lg', [1 1])
maximum= fscanf(fid, '%lg', [1 1])
rightascension_of_source= fscanf(fid, '%lg', [1 1])
declination_of_source= fscanf(fid, '%lg', [1 1])
RAsource=fscanf(fid, '%lg', [1 1])
Decsource=fscanf(fid, '%lg', [1 1])
TRI = delaunay(RAsourcevector,Decsourcevector);
nface=size(TRI);
for i=1:nface(1) 
color(i)=(1/3*(powermatrix3(TRI(i,1))+powermatrix3(TRI(i,2))+powermatrix3(TRI(i,3))));
color(i)==abs(color(i));
end
f3=figure(3)
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on','ParallelLabel','on','LabelFormat','none','MeridianLabel','on','MLabelParallel','equator');patchm(Decsourcevector,RAsourcevector,'Faces',TRI,'FaceColor','flat','LineStyle','none','FaceVertexCData',color')
xlabel('Right Ascension');
ylabel('Declination');
title({'SNR OVER THE SKY',' freq=',frequency3,'Hz,  ', '   mean power of noise =',meanofnoise,'   SNR =', SNR, '   stdv of noise = ', stdofnoise,'maxpower =', maximum,'RA,DEC(deg) =', rightascension_of_source, declination_of_source  });
colorbar('location','southoutside');
plotm([declination_of_source],[rightascension_of_source],'x');
plotm([Decsource*360/2/pi],[RAsource*360/2/pi],'+');
fclose(fid);
saveas(f3,'f3.fig')
