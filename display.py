from pylab import *

posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))
fake_data = loadtxt('fake_data.txt')

num = 200
L = 2./num
x = linspace(-1. + 0.5*L, 1. - 0.5*L, num)
y = x.copy()
[x, y] = meshgrid(x, y)
y = y[::-1]
r = sqrt(x**2 + y**2)

ion()
hold(False)

for i in xrange(0, posterior_sample.shape[0]):
	rho, rc, gamma = posterior_sample[i,0], posterior_sample[i,1], posterior_sample[i,2]


	f = rho/(1. + r/rc)**gamma
	subplot(1,3,1)
	imshow(fake_data)
	title('Data')
	gca().set_xticks([])
	gca().set_yticks([])

	subplot(1,3,2)
	imshow(f)
	gca().set_xticks([])
	gca().set_yticks([])

	title('Model %g'%(i+1))
	subplot(1,3,3)
	imshow(fake_data - f)
	title('Residuals')
	gca().set_xticks([])
	gca().set_yticks([])

	draw()

ioff()
show()

