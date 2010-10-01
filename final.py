"""Implementing the FrFT function """


import numpy as np
import scipy as sp
import matplotlib as mt
import matplotlib.pylab as pl
from pylab import zeros
from scipy import array





def bizint(y):
	N = len(y)
	y=np.mat(y)
	im=0

	if sp.sum(sp.absolute(sp.imag(y)))>0:
		im=1
		#print im
		imx= sp.imag(y)
		#print imx
		y=sp.real(y)
		#print y
	x2=np.array([[y],[sp.zeros((1,N))]])
	x2 = x2.T
	x2=x2.flatten()
	
	xf=sp.fft(x2)	

	if N%2 == 1:
		N1 = np.fix(N/2+1)
		N2 = 2*N-np.fix(N/2)
		#print N1,N2
		a = xf[0:N1]
		b = xf[N2:2*N]
		m = sp.zeros(N)
		g1=np.append(a,m)
		g = np.append(g1,b)
		#print g
		#print len(g)
		xint = 2*sp.real(sp.ifft(g))
		#print xint
		#print len(xint)	
		#pl.plot(xint)
		#pl.show()	
	else:
		
		a=xf[0:N/2]
		b=xf[2*N-N/2:2*N]
		print len(a),len(b),len(sp.zeros((N,1)))
		m = sp.zeros(N)
		g1 = np.append(a,m)
		#print f
		g = np.append(g1,b)
		#print g
		xint= 2*sp.real(sp.ifft(g))		
		#print xint
		
	if im==1:
		x2 = np.mat(imx)
		x2 = np.array([[x2],[sp.zeros((1,N))]])		
		x2 = x2.T
		x2 = x2.flatten()
		xf = fft(x2)
		if N%2 == 1:
			N1 = np.fix(N/2+1)
			N2 = 2*N-np.fix(N/2)
			c = xf[0:N1]
			d = xf[N2:2*N]
			e = sp.zeros(N)
			h1 = np.append(c,e)
			h = np.append(h1,d)
			xmint = 2*sp.real(sp.ifft(h))
		else:
			c = xf[0:N/2]
			d = xf[2*N-N/2:2*N]
			e = sp.zeros(N)
			h1 = np.append(c,e)
			h = np.append(h1,d)
			xmint = 2*sp.real(sp.ifft(h))
		xint = xint + 1j*xmint	


	xint = xint.flatten()
	return xint 

def chirpmethod(fc,a):
	deltax = np.sqrt(len(fc))
	#print deltax
	phi = a* np.pi/2
	N = np.fix(len(fc))
	alpha = 1/np.tan(phi)
	beta = 1/np.sin(phi)
	
	x = np.arange(-1*np.ceil(N/2),np.fix(N/2))/deltax
	#print x
	fc = fc.flatten()
	fc = fc[0:N]
	f1 = np.exp(-1j*np.pi*np.tan(phi/2)*x*x)  #Multiplying with chirp function"
	#print f1
	#pl.plot(f1)
	#pl.show()
	f1 = f1.flatten()
	fc = fc*f1
	#print fc
	#pl.plot(fc)
	#pl.show()
	x = x.flatten()
	
	t = np.arange(-N,N)/deltax
	hlptc = np.exp(1j*np.pi*beta*t*t)
	
	hlptc = hlptc.flatten()
	
	N2 = len(hlptc)
	#print N2
	p = np.ceil(np.log(N2+N-1)/np.log(2))
	#print p
	N3 = 2**p
	#print N3

	#print hlptc
	#print sp.zeros((N3-N,1))	
	#print len(fc)
	#print N2
	fc = np.append(fc,sp.zeros(len(fc)))
	
	hlptcz = np.array([hlptc,sp.zeros(N2)])
	#print len(hlptcz)
	
	fcz = np.array([fc,sp.zeros(len(fc))])
	#print len(fcz)

	#fcz1 = np.append(fcz,sp.zeros(len(fc)))	
	#print len(fcz1)

	
	Hcfft = sp.ifft(sp.fft(fcz)*sp.fft(hlptcz)) #convolution with chirp
	Hcfft = Hcfft.flatten()
	#print N, len(Hcfft)

	Hc = Hcfft[N:2*N]
	#print Hc
	#Hc = Hc.flatten()
	Aphi = np.exp(-1j*(np.pi*sp.sign(np.sin(phi))/4-phi/2))/np.sqrt(sp.absolute(np.sin(phi)))

	xx = np.arange(-np.ceil(N/2),np.fix(N/2)-1)/deltax
	f1 = f1.flatten()
	#print Aphi
	f1 = f1*Aphi
	#print len(f1),len(Hc)
	res = (f1*Hc)/deltax
	
	if np.fix(N/2)!=N/2:
		res2[1:N-1] = res[2:N]
		res2[-1:] = res[1]
		res = res2
	res = res.flatten()
	#pl.plot(res)
	#pl.show()
	#print res
	return res
	



def frft(y,a):

	N=len(y)
	y=y.flatten(1)
	a=1

	fc=bizint(y)
	#print fc
	#N1=len(fc)
	#print N1
	z = zeros(N)
	fc = np.append(z,fc)
	fc = np.append(fc,z)
	print fc
	#print fc
	print len(fc)
	#pl.plot(fc)
	#pl.show()

	flag = 0

	if a>0 and a<0.5:
		flag = 1
		a = a-1

	if a>-0.5 and a<0:
		flag = 2
		a = a+1

	if a>1.5 and a<2:
		flag = 3
		a = a+1

	if a>-2 and a<-1.5:
		flag = 4
		a = a+1

	res = fc

	if flag == 1 or flag == 3:
		res = chirpmethod(fc,1)

	if flag == 2 or flag == 4:
		res = chirpmethod(fc,-1)

	if a == 0:
		res = fc

	elif a == 2 or a == -2:
		res = fc[::-1,:]

	else:
		res = chirpmethod(res,a)


	res = res[N:3*N]
	#print len(res)
	#s = np.arange(1,len(res),2)
	#pl.plot(res)
	#pl.show()
	res = res[np.arange(1,len(res))]
	print res
	res = res.flatten()
	pl.plot(res)
	pl.show()
	return res

def main():
	a = 0.7
	t = np.arange(0,1.023,.001)
	f1 = 99
	f2 = 100
	f3 = 101

	sig1 = np.sin(2*np.pi*f1*t)
	sig2 = np.sin(2*np.pi*f2*t)
	sig3 = np.sin(2*np.pi*f3*t)

	sig = sig1+sig2+sig3 # blind source creation
	
	sigFr = frft(sig,a)
	
	#pl.plot(sigFr)
	#pl.show()
	fr1 = sp.real(sigFr)/np.argmax(np.real(sigFr))
	fr2 = sp.imag(sigFr)/np.argmax(np.imag(sigFr))
	
	#pl.plot(fr2)
	#pl.show()

if __name__ == '__main__':
	main()
