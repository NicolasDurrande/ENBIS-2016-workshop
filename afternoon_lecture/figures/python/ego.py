import numpy as np
import pylab as pb
from scipy.stats import norm
import GPy
from matplotlib import gridspec

pb.ion()

pb.rc('font', family='serif', size=15)
pb.rcParams['xtick.major.pad']='12'
pb.rcParams['ytick.major.pad']='12'
pb.rcParams['figure.figsize'] = 9, 5

def ftest(x):
	return(-np.sin(x*12)+x*1.2+2)

def ftestn(x):
	return(-np.sin(x*12)+x*1.2+2+np.random.normal(0,.1,x.shape))

x = np.linspace(0,1,101)[:,None]

#####################################################
# data
X = np.array([[0,.25,.5,.66,.9]]).T
X = np.array([[0,.3,.6,.9]]).T
F = ftest(X)
d = X.shape[1]

fig, ax = pb.subplots()
pb.plot(X,F,'kx',mew=2.5)
pb.xlim((-0.1,1.1))
pb.plot([0,1],[np.min(F)]*2,'k--',linewidth=1.5)
fig.subplots_adjust(bottom=0.15)
#fig.subplots_adjust(bottom=0.15)
pb.savefig('ego_obs.pdf',bbox_inches='tight')

#####################################################
#model 
# define noise variance
tau2 = 0

# define a kernel
kern = GPy.kern.Matern52(input_dim=d,variance=10,lengthscale=.4)

# define a model
m = GPy.models.gp_regression.GPRegression(X, F, kern)
m['.*noise'].fix(tau2)		# fix the noise variance
print m

m.plot(plot_limits=[0,1])
pb.plot([0,1],[np.min(F)]*2,'k--',linewidth=1.5)
pb.savefig('ego_modelInit.pdf',bbox_inches='tight')

#####################################################
# EI
def EI(x,m):
	mean, var = m.predict(x)
	var[var<0] = 0
	u = (np.min(m.Y) - mean)/(np.sqrt(var))
	ei = np.sqrt(var) * (u * norm.cdf(u) + norm.pdf(u))
	ei[np.isnan(ei)] = 0
	return(ei)

m = GPy.models.gp_regression.GPRegression(X, F, kern)
m['.*nois'] = 0
ei = EI(x,m)
xstar = x[np.argmax(ei)]
pb.figure(figsize=(9,7))
# plot model
gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1],width_ratios=[.25,20],) 
ax0 = pb.subplot(gs[1])
m.plot(plot_limits=[0,1],ax=ax0)
ax0.plot([0,1],[np.min(F)]*2,'k--',linewidth=1.5)
#ax0.plot(xstar,ftest(xstar),'rx',mew=2.5)
ax0.set_xticklabels([])
ax0.set_ylabel('$model$', fontsize=25)
ax0.get_yaxis().set_label_coords(-0.13,0.5)
ax0.legend_.remove()

# plot EI
ax1 = pb.subplot(gs[3])
ax1.plot(x,ei,'r',linewidth=2)
ax1.plot([xstar]*2,[0,np.max(ei)],'r--',linewidth=1.5)
ax1.set_ylim([-0.02, .25])
# mylocs = ax0.xaxis.get_majorticklocs()+[xstar]
# mylabs = ax0.xaxis.get_majorticklabels() + ['x']
# ax1.set_xticks(mylocs, mylabs)
ax1.locator_params(axis='y',nbins=2)
ax1.set_ylabel("$expected$ \n $improvement$", fontsize=25)
ax1.get_yaxis().set_label_coords(-0.1,0.5)
pb.savefig('ego_EI0.pdf',bbox_inches='tight')
#pb.savefig('ego_EI0bis.pdf',bbox_inches='tight')

# update DoE
X = np.vstack((X,xstar))
F = np.vstack((F,ftest(xstar)))

# new_label = '{0} xx'.format(xstar)
# ax1.set_xticks(xstar,new_label)

for i in range(1,6):
	m = GPy.models.gp_regression.GPRegression(X, F, kern)
	m['.*nois'] = 0
	ei = EI(x,m)
	xstar = x[np.argmax(ei)]
	pb.figure(figsize=(9,7))
	# plot model
	gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1],width_ratios=[.25,20],) 
	ax0 = pb.subplot(gs[1])
	m.plot(plot_limits=[0,1],ax=ax0)
	ax0.plot([0,1],[np.min(F)]*2,'k--',linewidth=1.5)
	#ax0.plot(xstar,ftest(xstar),'rx',mew=2.5)
	ax0.set_xticklabels([])
	ax0.set_ylabel('$model$', fontsize=25)
	ax0.get_yaxis().set_label_coords(-0.13,0.5)
	ax0.legend_.remove()
	# plot EI
	ax1 = pb.subplot(gs[3])
	ax1.plot(x,ei,'r',linewidth=2)
	ax1.plot([xstar]*2,[0,np.max(ei)],'r--',linewidth=1.5)
	ax1.set_ylim([-0.02, .25])
	ax1.locator_params(axis='y',nbins=2)
	#ax1.set_ylabel("$am\\'elioration$ \n $esp\\'er\\'ee$", fontsize=25)
	ax1.set_ylabel("$expected$ \n $improvement$", fontsize=25)
	ax1.get_yaxis().set_label_coords(-0.1,0.5)
	#save plot
	pb.savefig('ego_EI%i.pdf'%i,bbox_inches='tight')
	# update DoE
	X = np.vstack((X,xstar))
	F = np.vstack((F,ftest(xstar)))




K = m.kern.K(X,X)
lamb = np.linalg.eigvals(K)
round(lamb,2)
np.linalg.cond(K)


##############
## osborn taylor

X = np.linspace(0.5,0.51,2)[:,None]
K = m.kern.K(X,X)
np.linalg.cond(K)
np.linalg.eigvals(K)

Kd = K.copy()
Kd[0,1] = (K[0,1]-K[0,0])/.01
Kd[1,0] = (K[0,1]-K[0,0])/.01
Kd[1,1] = (K[1,1]+K[0,0]-2*K[0,1])/.0001
np.linalg.cond(Kd)
np.linalg.eigvals(Kd)

pb.figure(figsize=(5,5))
pb.plot(X,ftest(X),'kx',mew=1.5)
pb.xlim((0,1))
pb.ylim((0,6))
pb.savefig('osborn0',bbox_inches='tight')

pb.figure(figsize=(5,5))
pb.plot(X[0],ftest(X[0]),'kx',mew=1.5)
p = (ftest(X[1]) - ftest(X[0]))/.01
D = .05
pb.plot([X[0]-D,X[0]+D],[ftest(X[0])-p*D,ftest(X[0])+p*D],'r',linewidth=1.5)
pb.xlim((0,1))
pb.ylim((0,6))
pb.savefig('osborn1',bbox_inches='tight')

#####################################################
# noisy EI

def EIn1(x,m):
	mX, vX = m.predict(m.X)
	mn = GPy.models.gp_regression.GPRegression(X, mX, kern)
	mn['.*nois'] = 0.0
	mean, var = mn.predict(x)
	var[var<0] = np.inf
	u = (np.min(mX) - mean)/(np.sqrt(var))
	ei = np.sqrt(var) * (u * norm.cdf(u) + norm.pdf(u))
	return(ei)

def EMI(x,m):
	mean, var = m.predict(x)
	mX, vX = m.predict(m.X)
	var[var<0] = np.inf
	u = (np.min(mX) - mean)/np.sqrt(var+m['.*noise'])
	sig2 = var/(var+m['.*noise'])
	ein = np.sqrt(var+m['.*noise']) * (u*norm.cdf(u) + sig2*norm.pdf(u))
	return(ein)


x = np.linspace(0,1,100)[:,None]  
ei = EMI(x,m)

## llop
x = np.linspace(0,1,500)[:,None]  
X = np.linspace(0,.9,4)[:,None]
F = ftestn(X)

for i in range(6):
	m = GPy.models.gp_regression.GPRegression(X, F, kern)
	m['.*nois'] = 0.01
	ei = EIn1(x,m)
	xstar = x[np.argmax(ei)]
	mX, vX = m.predict(m.X)
	mn = GPy.models.gp_regression.GPRegression(X, mX, kern)
	mn['.*nois'] = 0.0
	pb.figure(figsize=(8,5))
	ax=pb.subplot(111)
	m.plot(plot_limits=[0,1],ax=ax)
	mn.plot(plot_limits=[0,1],ax=ax,linecol='g', fillcol='g')
	pb.plot(x,ei*20,'r',linewidth=2)
	pb.ylim((-.1,6))
	X = np.vstack((X,xstar))
	F = np.vstack((F,ftestn(xstar)))
	pb.plot([xstar]*2,[0,max(ei)*20],'r--',linewidth=1.5)
	pb.savefig('ego_EI1n%i.pdf'%i,bbox_inches='tight')


for i in range(12):
	m = GPy.models.gp_regression.GPRegression(X, F, kern)
	m['.*nois'] = 0.01
	ei = EMI(x,m)
	xstar = x[np.argmax(ei)]
	mX, vX = m.predict(m.X)
	pb.figure(figsize=(8,5))
	ax=pb.subplot(111)
	m.plot(plot_limits=[0,1],ax=ax)
	pb.plot(x,ei*20,'r',linewidth=2)
	pb.ylim((-.1,6))
	X = np.vstack((X,xstar))
	F = np.vstack((F,ftestn(xstar)))
	pb.plot([xstar]*2,[0,max(ei)*20],'r--',linewidth=1.5)
	# pb.savefig('ego_EI1n%i.pdf'%i,bbox_inches='tight')


m.predict(np.array([[.13]]))

ones = np.ones((20,20))

ones - np.dot(np.dot(ones,np.linalg.inv(ones+np.eye(20))),ones)

#####################################################
# inversion

# load data
X = np.array([[0,.25,.5,.66,.9]]).T
F = ftest(X)
pb.plot(X,F,'kx')
d = X.shape[1]

# define noise variance
tau2 = 0

# define a kernel
kern = GPy.kern.Matern52(input_dim=d,variance=10,lengthscale=.4)

# define a model
m = GPy.models.gp_regression.GPRegression(X, F, kern)
m['.*noise'].fix(tau2)		# fix the noise variance

m.plot(plot_limits=[0,1])
pb.plot([0,1],[3.2]*2,'k--',linewidth=1.5)
pb.savefig('inv.pdf',bbox_inches='tight')


##
mean, var = m.predict(x)
prob = norm.pdf(3.2,mean,var)
m.plot(plot_limits=[0,1])
pb.plot([0,1],[3.2]*2,'k--',linewidth=1.5)
pb.plot(x,prob)
pb.ylim((-.1,6))
pb.savefig('invproba.pdf',bbox_inches='tight')


xstar = x[np.argmax(prob)]
X = np.vstack((X,xstar))
F = np.vstack((F,ftestn(xstar)))
m = GPy.models.gp_regression.GPRegression(X, F, kern)
m['.*noise'].fix(tau2)		# fix the noise variance

mean, var = m.predict(x)
prob = norm.pdf(3.2,mean,var)
m.plot(plot_limits=[0,1])
pb.plot([0,1],[3.2]*2,'k--',linewidth=1.5)
pb.plot(x,prob)
pb.ylim((-.1,6))
pb.savefig('invproba1.pdf',bbox_inches='tight')


xstar = x[np.argmax(prob[:250])]
X = np.vstack((X,xstar))
F = np.vstack((F,ftestn(xstar)))
m = GPy.models.gp_regression.GPRegression(X, F, kern)
m['.*noise'].fix(tau2)		# fix the noise variance

mean, var = m.predict(x)
prob = norm.pdf(3.2,mean,var)
m.plot(plot_limits=[0,1])
pb.plot([0,1],[3.2]*2,'k--',linewidth=1.5)
pb.plot(x,prob)
pb.ylim((-.1,6))
pb.savefig('invproba2.pdf',bbox_inches='tight')


xstar = x[np.argmax(prob[:250])]
X = np.vstack((X,xstar))
F = np.vstack((F,ftestn(xstar)))
m = GPy.models.gp_regression.GPRegression(X, F, kern)
m['.*noise'].fix(tau2)		# fix the noise variance

mean, var = m.predict(x)
prob = norm.pdf(3.2,mean,var)
m.plot(plot_limits=[0,1])
pb.plot([0,1],[3.2]*2,'k--',linewidth=1.5)
pb.plot(x,prob)
pb.ylim((-.1,6))
pb.savefig('invproba3.pdf',bbox_inches='tight')


