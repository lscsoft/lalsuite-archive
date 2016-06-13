#
#  Compute Feldman-Cousins upper limit values for
#  Gaussian distribution and save in file fc.txt
#  
#  You need quanta to run this
#
#  See: Gary J. Feldman, Robert D.Cousins
#  A Unified Approach to the Classical Statistical Analysis of Small Signals
#  arXiv:physics/9711021
#

global alpha;

#
# Confidence level

# 95% Feldman-Cousins limits
alpha=0.95;

# 90% Feldman-Cousins limits
#alpha=0.9;

# Range and step of mu values to sample
global mu;
mu=0.0;
mu_max=9.0;
mu_step=0.0005;

#
# Range and step of x values to sample
#
x_min=-5;
x_step=0.001;
x_max=5;

#
# The R function as in FC paper
#

function y=R(x, mu)
 if(x>=0)y=exp(-(x-mu)**2/2);
 	else y=exp(x*mu-mu**2/2);
 end
end

#
#  Normal distribution function with mean mu
#
function y=NDIST(x, mu)
y=(erf((x-mu)/sqrt(2.0))+1.0)/2.0;
end

#
# Inverse of normal distribution function with mean mu,
# i.e. quantiles
#
function y=NDISTINV(x, mu)
y=erfinv(2*x-1)*sqrt(2.0)+mu;
end


#
# Right Bound function
# Given mu and x1 compute x2 from condition
# NDIST(x2, mu)-NDIST(x1,mu)=alpha
#
function y=RB(x, mu, alpha)
a=NDIST(x, mu)+alpha;
if(a>1.0)y=10^6;
	else
	y=NDISTINV(a, mu);
	end;
end

#
# This function is the one actually passed to fsolve
#
function y=f(x)
global alpha
global mu
y=R(x, mu)-R(RB(x, mu, alpha), mu);
end

k=1;

while(mu<mu_max)

if(abs(mu-0.6)>0.02)
	[x, info]=fsolve("f",mu-0.6);
	else
	[x, info]=fsolve("f", mu-0.4);
	endif

MU(k)=mu;
UL(k)=x;
LL(k)=RB(x, mu, alpha);
k++;
mu=mu+mu_step;
endwhile

#plot(MU,UL)
#plot(MU,LL)
#plot([LL',UL'],MU)


#
#
# Invert MU->[LL', UL'] mapping to obtain actual confidence intervals
#

x=x_min;
k=1;
while(x<x_max)
X(k)=x;
[val, key]=min(UL<x);
MU_UL(k)=MU(key);
[val, key]=min(LL<x);
MU_LL(k)=MU(key);
x=x+x_step;
k++;
end

plot(X, [MU_LL', MU_UL'])

#
# Save limits in the file
#
LIMITS=[X', MU_LL', MU_UL'];
save("fc.txt", "LIMITS")

