from pylab import *

seed(2)

num = 200
L = 2./num
x = linspace(-1. + 0.5*L, 1. - 0.5*L, num)
y = x.copy()
[x, y] = meshgrid(x, y)
y = y[::-1]
r = sqrt(x**2 + y**2)

F = 1.
rc = 0.3
gamma = 2.

f = 1. + F/(2.*pi*rc**2)*exp(-0.5*(r/rc)**2)
substructure = zeros(f.shape)


for i in xrange(0, 1000):
	xc = -1. + 2.*rand()
	yc = -1. + 2.*rand()
	rr = (x - xc)**2 + (y - yc)**2
	substructure += (-0.2*log(rand()))*exp(-0.5*rr/0.02**2)

f += substructure
f += sqrt(0.1*f)*randn(num, num)
savetxt('fake_data.txt', f)

subplot(1,2,1)
imshow(f)

subplot(1,2,2)
imshow(substructure)
show()

