from pylab import *

posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('data.txt')
valid = data > -1E250

num = 200
L = 2./num
x = linspace(-1. + 0.5*L, 1. - 0.5*L, num)
y = x.copy()
[x, y] = meshgrid(x, y)
y = y[::-1]

ion()
hold(False)

for i in xrange(0, posterior_sample.shape[0]):
	rho, rc, gamma, xc, yc, q, theta = posterior_sample[i,0], posterior_sample[i,1],\
			 posterior_sample[i,2], posterior_sample[i, 3], posterior_sample[i, 4]\
			, posterior_sample[i, 5], posterior_sample[i, 6]

	xx =  cos(theta)*(x - xc) + sin(theta)*(y - yc);
	yy = -sin(theta)*(x - xc) + cos(theta)*(y - yc);
	r = sqrt(q*xx**2 + yy**2/q);
	f = rho*exp(-(r/rc)**gamma)
	subplot(1,3,1)
	imshow(data*valid)
	title('Data')
	gca().set_xticks([])
	gca().set_yticks([])

	subplot(1,3,2)
	imshow(f)
	gca().set_xticks([])
	gca().set_yticks([])

	title('Model %g'%(i+1))
	subplot(1,3,3)
	imshow(data*valid - f*valid)
	title('Residuals')
	gca().set_xticks([])
	gca().set_yticks([])

	draw()

ioff()
show()


