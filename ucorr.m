%ucorr
%Incertezza variabili correlate
%Zhang 2006 Metrologia 
pkg load signal
pkg load io

close all;
clear all;



%x = randn(200,1);
%for count = 1:100, x(2*count) = x(2*count-1); end;
%x = cumsum(x);
%x = x + 10; %mean
%x = input('x:');
%x = load('-ascii','6010B_10_1.dat');

%col = 'B'
%rowstart = 84;
%rowend = 333; % 819;
%cut = 20;
%xfull = xlsread ('QHR_4-10-1k_STD_20221110T1210_50uA.xlsx', 'QHR_4-10-1k_STD_20221110T1210_5', sprintf('%c%d:%c%d',col,rowstart,col,rowend));
%x = xfull(cut+1:length(xfull));

cut = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Read .MEA files

[fname, fpath, fltidx] = uigetfile ();
fnamecomp = [fpath fname];
fileid = fopen(fnamecomp);

%hdrRows = 83;
%hdrData = textscan(fileid,'%s',hdrRows, 'Delimiter','\n');
%matData = textscan(fileid,'%f %f %s', 'Delimiter',';','CollectOutput',true);

lines = textscan(fileid,'%s','delimiter','\n'); %this is read as a cell array of cell array of strings

lines=lines{1}; %reduce to a single cell array of strings

%run through all lines to find two instances of ';'
countgoodlines = 0;
ratio = [];
for linek = 1:length(lines),
  if (sum(char(lines(linek)) == ';') >= 2), % 6010D has 2 ";" with no temperature and 3 with; the 6020Q has 3 ";" 
      data = sscanf(char(lines(linek)),'%f;%f'); %Read the first two numbers in the rows
      if length(data) > 1, %Good lines
        if (mod(data(1),1) == 0) ratiok = data(2); %if it is an integer, the file includes the number of sample column
          else ratiok = data(1);
        end;
        countgoodlines=countgoodlines+1;
        ratio(countgoodlines) = ratiok;
      end;
   end;   
end;
    
fclose(fileid);

xfull = ratio';

%%%%%%%%%%%%%%%%%%%%%%%%



x = xfull(cut+1:length(xfull));

N = length(x);

meanx = mean(x);
stdx = std(x);
stdmeanx = std(x)./sqrt(N);

[R,lag] = xcov(x);
%shorten the vector to lag >0
R0 = R(N); %Special case lag 0
R = R((N+1):(2*N-1)); %From lag 1

rho = R./R0; %normalized

rhoband = zeros(N-1,1);
rhoband(1) = 1.96*sqrt(1/N);
for i=2:N-1, 
  rhoband(i) = 1.96*sqrt( (1 + 2*sum(rho(1:(i-1)).^2))./N );
  %rhoband(i) = 3*sqrt( (1 + 2*sum(rho(1:(i-1)).^2))./N );
end;

Pc = find((abs(rho)>rhoband));

Nc = max(find(abs(rho)>rhoband));
Nr = min(Nc, floor(N/4));

ux_nocorr = std(x)./sqrt(N);
srho = 0;
%for i = 1:Nr,
for i = Pc',
  srho = srho + (N-i).*rho(i);
  %disp(sprintf('\n %d %f %f ',i,rho(i), srho));
end;

ux_expfactor = sqrt( 1 + (2/N).*srho);
ux_corr = ux_expfactor * ux_nocorr;

figure(1);
plot(-(cut-1):N,xfull-meanx,'c',1:N,x-meanx,'b','linewidth',3);
set(gca, "fontsize", 14);
grid on;
zoom on;


figure(2)
hold on;
plot(1:(N-1),rhoband,'c-','linewidth', 3);
plot(1:(N-1),-rhoband,'c-', "linewidth", 3);
plot(1:(N-1),rho,'bo')
plot(Nc+1,rho(Nc+1),'r*','MarkerSize',20);
hold off;
set(gca, "fontsize", 14);
grid on;
zoom on;


disp(sprintf('mean(x)            =%14.11f', meanx));
disp(sprintf('std(x)             =%e', stdx));
disp(sprintf('stdmean(x)         =%e', stdmeanx));
disp(sprintf('stdmeanrel(x)      =%e', stdmeanx./meanx));
disp(sprintf('ux(x)              =%e', ux_corr));
disp(sprintf('uxrel(x)           =%e', ux_corr./meanx));
disp(sprintf('correlation length =%e', max(Pc)));
disp(sprintf('corr exp factor    =%f', ux_expfactor));
disp(sprintf(''));
disp(sprintf('%s,%12.9f,%3.2e',fname,meanx,ux_corr./meanx));
