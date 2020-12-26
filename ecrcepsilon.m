function G=ecrcepsilon(Gij)
global r epsilon_a epsilon_b f error
x=r*norm(Gij);
if norm(Gij)<error
    G=1/epsilon_a*f+1/epsilon_b*(1-f);
else
    G=(1/epsilon_a-1/epsilon_b)*f*2*besselj(1,x)/x;% expansion coefficient of epsilon in reciprocal space
end
end