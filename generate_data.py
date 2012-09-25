from pylab import *

seed(0)

num = 200
L = 2./num
x = linspace(-1. + 0.5*L, 1. - 0.5*L, num)
y = x.copy()
[x, y] = meshgrid(x, y)
y = y[::-1]
r = sqrt(x**2 + y**2)

rho = 1.
rc = 0.3
gamma = 2.

f = rho/(1. + r/rc)**gamma

for i in xrange(0, 1000):
	xc = -1. + 2.*rand()
	yc = -1. + 2.*rand()
	rr = (x - xc)**2 + (y - yc)**2
	f = f + 0.1*rand()*exp(-0.5*rr/0.02**2)


f = f + 0.03*randn(num, num)
savetxt('fake_data.txt', f)
imshow(f)
show()

