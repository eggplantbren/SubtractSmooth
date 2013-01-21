from pylab import *

posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))
data = loadtxt('fake_data.txt')
valid = data > -1E250

num = 200
L = 2./num
x = linspace(-1. + 0.5*L, 1. - 0.5*L, num)
y = x.copy()
[x, y] = meshgrid(x, y)
y = y[::-1, :]

ion()
hold(False)

for i in xrange(0, posterior_sample.shape[0]):
	f = posterior_sample[i, 13:].reshape((num, num))

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


