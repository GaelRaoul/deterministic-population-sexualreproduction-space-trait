clear all; %close all;

%octave --gui

Kc=20;
cREC=[];vS=[];vN=[];vZS=[];vZN=[];IkREC1=[];IkREC2=[];lambdaREC=[];AREC=[];BREC=[];CREC=[];

%fid = fopen('vRECsexueRonceVLE01.csv','wt');

#Iterations for different values of the climate change speed c
for kc=14:20

kc
  m=400;
  n=2^10;
  tfinal=101;
  L=350;
  kplot=round(tfinal);
  lplot=20;
  rhoin=.6;

  sigma=1;rmax=1;Vs=1;kmodel=1;eta=1;
  b=0.1
  VLE=0.3
  c=2.4*kc/Kc
  %c=1.5
  %limites VLE 0.1: 1.65
  %limites VLE 1: 5.5

filename3 = ["data/eachspeedVLEi" num2str(VLE) "b" num2str(b) "tfinal" num2str(tfinal)  "kc" num2str(kc)  "c" num2str(c) ".csv"];
fid3 = fopen(filename3,'w');
fprintf(fid3,'t,c,vS,vN,vZS,vZN,xS,xN,Ik1,Ik2,lambda,VLE,b,C,tfinal\n');

  alpha=VLE/(rmax*Vs); beta=sqrt(sigma^2/(2*VLE*rmax))*b;gamma=eta;nu=sqrt(2/(rmax*sigma^2))*c;

  %alpha=0.05;
  %beta=2.8289;
  %gamma=1;nu=0.5;

##  alpha=0.05;
##  beta=1.5;
##  gamma=1;nu=0;

  r_ast=1-VLE/(2*rmax*Vs)
  AKB=VLE/(r_ast*rmax*Vs)
  BKB=b*sigma/(2*r_ast*rmax*sqrt(Vs))
  climatespeed=sqrt(2/(rmax*r_ast*sigma^2))*c
  cUD=AKB*sqrt(2)/(BKB*sqrt(2/(rmax*r_ast*sigma^2)))
  %cUD2=VLE*(r_ast)^(3/2)*sqrt(rmax)/b
  cUD3=2*VLE*rmax*sqrt(Vs-VLE/(2*rmax^2))/b
  cPEase=2*sqrt(1-(BKB*sqrt(2)+AKB)/2)/sqrt(2/(rmax*r_ast*sigma^2))
  cPease2=sqrt(2)/b*sqrt((rmax-sqrt(b^2*sigma^2/Vs)/2)*b^2*sigma^2)

  K=L*beta+40;

  #Phase d'initialisation
  dx=L/(m); x=(1:m)*dx;
  dv=K/n; v=(1:n)*dv;
  %dt=0.2*dx^2/2;
  dt=0.1

  #Inintialisation de la population
  nd=zeros(m,n);
  Zopt0=K/2+beta*(x-L/2); Zopt=Zopt0;prop=0.5;
  for j=1:n
      for i=1:m
          nd(i,j)=(0.05/dv)*exp(-(abs(x(i)-L*prop))^2/20-(v(j)-Zopt(i))^2/2);%+(v(j)-Zopt(i))^2));
      end
  end
  nd0=nd; ndd=nd;
  Iplot=zeros(kplot,m);
  nplot=zeros(kplot,m,n);
  R=zeros(m,n);

  #Calculation of the reproduction kernel to compute the reproduction term as a double convoilution
  clear v2
  for i=1:n
      v2(i)=-v(n-i+1);
  end
  v2=[v2 0 v];
  Q=exp(-(v2).^2/4)/sqrt(pi);


% interior size in first direction
rdiff=dt/(dx*dx);
Mdiff = m - 2;
e = ones(Mdiff,1);
% matrix for implicit diffusion
A = spdiags([-rdiff*e, (1+2*rdiff)*e, -rdiff*e], [-1,0,1], Mdiff, Mdiff);


  #Initialisation des indices et tableaux
  xshift=0;
  yshift=0;
  xt=[];xtm=[];Zt=[];Ztm=[];Ikrec=[];
  tx=[];speedrec=[];speedrecm=[];speedrecz=[];speedreczm=[];maxrec=[];
  kpt=0;
  lpt=0;
  temps=0;

  #Boucle en temps
  while temps < tfinal%1.5*dt%

    Zopt=K/2+beta*(x+xshift-L/2-nu*temps);
    temps;
    dt=min(dt, tfinal-temps);
    nd00=nd;
    # Calculation of the time step iteration
    for i=1:m
     iopt=round((Zopt(i)-v(1))/(v(n)-v(1))*n);
     im=max(round(iopt-(10+sqrt(2/alpha))/dv),1);ip=min(round(iopt+(10+sqrt(2/alpha))/dv),n);
     di=ip-im+1;
     Ik(i)=norm(nd(i,:),1)*dv;
     clear tric truc chose;
     # calculation of the birth term
     if Ik(i)>10^-5
      f=nd(i,im:ip);
      Qtemp=Q(floor(length(Q)/2-di)+1:ceil(length(Q)/2+di));
      tric=conv(conv(Qtemp,f)*dv,f)*dv;
      for j=1:di truc(j)=tric(2*j+di-1);end
      chose=[nd(i,1:im-1) (gamma)/Ik(i)*truc nd(i,ip+1:n)];
     else
      chose=zeros(1,n);
     end
     ndd(i,:)=dt*chose+nd(i,:).*exp(-dt*((gamma-1)+min(alpha/2*(v-yshift-Zopt(i)).^2,20000)+Ik(i)));
    end
    #Dispersion
    %nd(2:m-1,2:n-1)=ndd(2:m-1,2:n-1)+(dt/(dx*dx))*(ndd(1:m-2,2:n-1)+ndd(3:m,2:n-1)-2*ndd(2:m-1,2:n-1));
    nd(2:m-1,2:n-1) = A \ ndd(2:m-1,2:n-1);
    #Boundary values
    nd(m,2:n-1)=nd(m-1,2:n-1);
    for j=2:n-10
     interp=dx*beta/dv;
     deltaj=floor(interp)+1;
     nd(1,j)=interp/deltaj*nd(2,j+deltaj)+(deltaj-interp)/deltaj*nd(2,j);
    end
    nd(1,n-10:n-1)=nd(2,n-10:n-1);
    nd(2:m-1,1)=0;
    nd(2:m-1,n)=0;
    nd(1,1)=.5*nd(2,1) +.5*nd(1,2);
    nd(1,n)=.5*nd(2,n) +.5*nd(1,n-1);
    nd(m,1)=0;
    nd(m,n)=.5*nd(m,n-1) +.5*nd(m-1,n);

    #Recording of the position of the population rightmost and leftmost)
    tx=[tx temps];Ikrec=[Ikrec norm(Ik,1)];

    #Calculation of the rightmost position occupied by the population
    k=m;thresholdI=0.03;
    while Ik(k)<thresholdI&&k>1
        k=k-1;
    end
    k=min(k, length(Ik)-1);
    thetaI=(thresholdI-Ik(k+1))/(Ik(k)-Ik(k+1));
    thetaI=min(1,thetaI);thetaI=max(0,thetaI);
    XtI=thetaI*k*dx+(1-thetaI)*(k+1)*dx+xshift;
    xt=[xt XtI];
    speed=(xt(length(tx))-xt(floor(length(tx)*0.9)+1))/(max(dt,tx(length(tx))-tx(floor(length(tx)*0.9)+1)));
    speedrec=[speedrec speed];
    zk=sum(nd(k,:).*v)/sum(nd(k,:));
    zk1=sum(nd(k+1,:).*v)/sum(nd(k+1,:));
    ZtI=thetaI*zk+(1-thetaI)*zk1-yshift;
    Zt=[Zt ZtI];
    speedz=(Zt(length(tx))-Zt(floor(length(tx)*0.9)+1))/(max(dt,tx(length(tx))-tx(floor(length(tx)*0.9)+1)));
    speedrecz=[speedrecz speedz];

    #Calculation of the leftmost position occupied by the population
    km=2;
    while Ik(km)<thresholdI&&km<m
        km=km+1;
    end
    thetaI=(thresholdI-Ik(km-1))/(Ik(km)-Ik(km-1));
    thetaI=max(0,thetaI);thetaI=min(1,thetaI);
    XtIm=thetaI*km*dx+(1-thetaI)*(km-1)*dx+xshift;
    xtm=[xtm XtIm];
    speedm=(xtm(length(tx))-xtm(floor(length(tx)*0.9)+1))/(max(dt,tx(length(tx))-tx(floor(length(tx)*0.9)+1)));
    speedrecm=[speedrecm speedm];
    zk=sum(nd(km,:).*v)/sum(nd(km,:));
    zk1=sum(nd(km-1,:).*v)/sum(nd(km-1,:));
    ZtIm=thetaI*zk+(1-thetaI)*zk1-yshift;
    Ztm=[Ztm ZtIm];
    speedzm=(Ztm(length(tx))-Ztm(floor(length(tx)*0.9)+1))/(max(dt,tx(length(tx))-tx(floor(length(tx)*0.9)+1)));
    speedreczm=[speedreczm speedzm];
  maxrec=[maxrec norm(Ik,Inf)];

    #Calculation populstion size evolution rate in log scale
    lambda=log(Ikrec(length(tx))/Ikrec(floor(length(tx)/2)+1))/(max(dt,tx(length(tx))-tx(floor(length(tx)/2)+1)));

    #After a sufficient number of step, the following instences are launched, to recenter the population and plot the results
    if (temps > tfinal*(kpt/kplot)|| temps>tfinal-2*dt)
      temps
      kpt=kpt+1;
      imin=m;imax=1;
      #Shift in x if the population approaches the large x limit of the calculation window
      Mx = sum(nd, 2);
      xcm=sum(x.*(Mx'))/sum(Mx);
      if xcm>0.6*x(length(x))
        nd1=nd;
        xshift=xshift+x(round(length(x)/5));
        k = round(length(x)/5);
        ndtemp=nd';
        ndtemp2 = [ndtemp(:,k+1:size(ndtemp,2)),zeros(size(ndtemp,1),k)];
        nd=ndtemp2';
        nd = [zeros(size(nd,1),k),nd(:,k+1:end)];
        tempsxshift=temps
        nd2=nd;
      end

      #Shift in y if the population approaches the large x limit of the calculation window
      My = sum(nd, 1);
      ycm=sum(v.*(My))/sum(My);
##      if ycm<0.4*v(length(v))
##        nd3=nd;
##        yshift=yshift+v(round(length(v)/5));
##        k = round(length(v)/5);
##        nd = [zeros(size(nd,1),k),nd(:,1:size(nd,2)-k)];
##        %nd = [zeros(k, size(nd,2)); nd(1:end-k,:)];
##        tempsyshift=temps
##        nd4=nd;
##      end

      #Calculation of the mean phenotypic trait
      for i=1:m
        if Ik(i)>0.001
          Zk(i)=norm(nd(i,1:n).*v,1)/norm(nd(i,1:n),1);
          imin=min(imin,i);imax=max(imax,i);
        else
          Zk(i)=0;
        end
      end
      #Graphical representations of the results
##          figure(1)
##          pcolor(x(10:m-10),v(10:n-10),nd(10:m-10,10:n-10)'),shading interp, drawnow;
##          figure(2)%plot(x,Zk,'b',x,Zopt,'r')
##          plot(x, Ik),title('N'), drawnow;
##          if imin<imax
##              figure(3)
##              plot(x(imin:imax),Zk(imin:imax),'r', x,Zopt,'k',x,Zopt0,'k--',x, Zopt-sqrt(2/alpha),'b',x, Zopt+sqrt(2/alpha)),title('Z'), drawnow;
##          end
##        figure(4)
##  	    plot(tx,xt,'k',tx,xtm,'r'),legend('northen tip','southern tip'),xlabel('t'),ylabel('x where I=0.5'),title(['Position of the front t->X(t), mean speed south=' num2str(speedm) 'north=' num2str(speed)]), drawnow;
##              figure(5)
##        semilogy(tx,Ikrec,'b',tx,exp(lambda*tx),'r' ),xlabel('t'),ylabel('total pop size')

      if temps > tfinal*(lpt/lplot)
        lpt=lpt+1;
      end
    end
    temps=temps+dt;
    fprintf(fid3,'%d,%d, %d, %d, %d, %d, %d, %d,%d,%d,%d,%d,%d,%d,%d\n',[temps nu*sqrt(sigma^2/2) sqrt(2/(sigma^2))*speedm sqrt(2/(sigma^2))*speed speedzm*sqrt(VLE*Vs) XtIm XtI speedz*sqrt(VLE*Vs) norm(Ik,Inf) norm(Ik,1) lambda AKB BKB climatespeed tfinal]');
    %fprintf(fid3,'%d\n',[temps ]');

end
fclose(fid3)
krec=k

  #Recording of the main sdata of the calculations
  cREC=[cREC nu*sqrt(sigma^2/2)];vS=[vS speedm*sqrt(sigma^2/2)];vN=[vN speed*sqrt(sigma^2/2)];vZS=[vZS speedzm*sqrt(VLE*Vs)];vZN=[vZN speedz*sqrt(VLE*Vs)];
  IkREC1=[IkREC1 Ikrec(floor(length(Ikrec)*0.8)+1)];IkREC2=[IkREC2 Ikrec(length(Ikrec))];
  lambdaREC=[lambdaREC lambda];AREC=[AREC alpha];BREC=[BREC beta];CREC=[CREC gamma];

end

filename2 = ["data/PropspeedsVLEi" num2str(VLE) "b" num2str(b) "tfinal" num2str(tfinal)  "kc" num2str(kc)  "c" num2str(c) ".csv"];
fid2 = fopen(filename2,'w');
fprintf(fid2,'c,vS,vN,vZS,vZN,Ik1,Ik2,lambda,VLE,b,C,tfinal\n');
for j=1:length(cREC)
  fprintf(fid2,'%d, %d, %d, %d, %d, %d, %d,%d,%d,%d,%d,%d\n',[cREC(j) vS(j) vN(j) vZS(j) vZN(j) IkREC1(j) IkREC2(j) lambdaREC(j) AREC(j) BREC(j) CREC(j) tfinal]');
end
fclose(fid2)


for i=1:m
  xk0(i)=x(i)*sqrt((sigma^2)/(2*rmax));
  xk(i)=(x(i)+xshift)*sqrt((sigma^2)/(2*rmax));
  Zoptk(i)=sqrt(VLE)*Zopt(i);
  Ik(i)=norm(nd(i,:),1)*dv;
  xshiftrec(i)=xshift;
  if Ik(i)>0.000001
    Zk(i)=sqrt(VLE)*norm(nd(i,1:n).*v,1)/norm(nd(i,1:n),1);
    Vk(i)=VLE*norm(nd(i,1:n).*(v-Zk(i)/sqrt(VLE)).^2,1)/norm(nd(i,1:n),1);
    imin=min(imin,i);imax=max(imax,i);
  else
    Zk(i)=0;
    Vk(i)=0;
  end
  Zoptk0(i)=sqrt(VLE)*Zopt0(i);
  Ik0(i)=norm(nd0(i,:),1)*dv;
  if Ik0(i)>0.000001
    Zk0(i)=sqrt(VLE)*norm(nd0(i,1:n).*v,1)/norm(nd0(i,1:n),1);
    Vk0(i)=VLE*norm(nd0(i,1:n).*(v-Zk0(i)/sqrt(VLE)).^2,1)/norm(nd0(i,1:n),1);
    imin=min(imin,i);imax=max(imax,i);
  else
    Zk0(i)=0;
    Vk0(i)=0;
  end
end

M = [xk(:), Ik(:),Zk(:),Vk(:),Zoptk(:), xk0(:), Ik0(:),Zk0(:),Vk0(:),Zoptk0(:),xshiftrec(:)];

filename = ["data/RECVLEi" num2str(VLE) "b" num2str(b) "tfinal" num2str(tfinal) "kc" num2str(kc) "c" num2str(c) ".csv"];
fid = fopen(filename,'w');
fprintf(fid,'x,I,Z,V,Zopt,x0,I0,Z0,V0,Zopt0,xshift\n');
fclose(fid);
dlmwrite(filename, M, '-append');









