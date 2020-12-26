% by Shuaifei Ren & Zhuoxiong Liu
clear
close all
global r error epsilon_a epsilon_b f
error=1.0e-6;
a=1;    %characteristic length of the lattice
c=1;    %velocity of light
a1=a*[1 0];a2=a*[0 1];  %the primitive lattice vectors of the direct lattice
b1=2*pi/a*[a2(2) -a2(1)];b2=2*pi/a*[-a1(2) a1(1)];  % the primitive lattice vectors of the reciprocal lattice
r=0.2*a;              % radius of the dielectric columns
f=pi*r^2/a^2;    %the filling fraction
epsilon_a=8.9;
epsilon_b=1.0;
Gamma=[0 0];Chi=[pi/a 0];M=[pi/a pi/a];     %Points of high-symmetry on the Brillouin zone
mode=0;     %0 for TE mode and 1 for TM mode

Nk=20;
GammaChi=zeros(Nk+1,2);
ChiM=zeros(Nk+1,2);
MGamma=zeros(Nk+1,2);
for i=1:(Nk+1)
    GammaChi(i,:)=Gamma+(i-1)*(Chi-Gamma)/Nk;
    ChiM(i,:)=Chi+(i-1)*(M-Chi)/Nk;
    MGamma(i,:)=M-(i-1)*(M-Gamma)/Nk;
end
kkvector={GammaChi,ChiM,MGamma};

N=10;
NG=(2*N+1)^2;
G=zeros(NG,2);
i=1;
for h1=-N:N
    for h2=-N:N
        G(i,:)=h1*b1+h2*b2;%reciprocal vectors
        i=i+1;
    end
end;

omega=zeros(NG,3*(Nk+1));
x=linspace(0,1,3*(Nk+1));

for t=1:3
    F=zeros(NG,NG);
    kvector=kkvector{t};
    for k=1:Nk+1
        for i=1:NG
            for j=1:NG
                if mode==0
                    F(i,j)=dot((kvector(k,:)+G(i,:)),(kvector(k,:)+G(j,:)))*ecrcepsilon((G(i,:)-G(j,:)));
                elseif mode==1
                    F(i,j)=norm(kvector(k,:)+G(i,:))*norm(kvector(k,:)+G(j,:))*ecrcepsilon((G(i,:)-G(j,:)));
                else
                    disp('ERROR! no such mode,mode should be 0 for TE mode or 1for TM mode')
                    return
                end
            end
        end
        omega(:,(t-1)*(Nk+1)+k)=sqrt(eig(F))*a/(2*pi);
    end
end

if mode==0
    for i=1:NG
        plot(x,omega(i,:),'r','Linewidth',1)
        hold on
    end
    text(5/63,0.35,'TE modes','color',[1 0 0]);
elseif mode==1
    for i=1:NG
        plot(x,omega(i,:),'color',[0.09 0.32 0.59],'Linewidth',1)
        hold on
    end
    text(5/63,0.08,'TM modes','color',[0.09 0.32 0.59]);
end

ylim([0 0.8])
plot([1/3 1/3], get(gca, 'YLim'),'k');
plot([2/3 2/3], get(gca, 'YLim'),'k');
set(gca,'xtick',[]);
ylabel('Frequency   \omegaa/2\pic');

text(0-0.018,-0.03,'\Gamma','FontSize',15);
text(1/3-0.018,-0.03,'\chi','FontSize',15);
text(2/3-0.018,-0.03,'M','FontSize',15);
text(1-0.018,-0.03,'\Gamma','FontSize',15);
