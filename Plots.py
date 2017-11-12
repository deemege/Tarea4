import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import write

c= np.genfromtxt("datosCuerda.txt").T
cP= np.genfromtxt("datosCuerdaPerturbada.txt").T
s= np.genfromtxt("sonido.txt").T

x,y0,y1,y2,y3=c[0], c[1], c[2], c[3], c[4]
plt.plot(x,y0,label="t=0")
plt.plot(x,y1,label="t=T/8")
plt.plot(x,y2,label="t=T/4")
plt.plot(x,y3,label="t=T/2")
plt.title("Cuerdas con extremos fijos")
plt.xlabel("x (m)")
plt.ylabel("Amplitud (m)")
plt.legend()
plt.savefig("cuerda.png")
plt.close()

y0P,y1P,y2P,y3P= cP[1], cP[2], cP[3], cP[4]
plt.plot(x,y0P,label="t=0")
plt.plot(x,y1P,label="t=T/8")
plt.plot(x,y2P,label="t=T/4")
plt.plot(x,y3P,label="t=T/2")
plt.title("Cuerda con un extremo fijo y el otro extremo perturbado")
plt.xlabel("x (m)")
plt.ylabel("Amplitud (m)")
plt.legend()
plt.savefig("cuerdaPerturbada.png")

audio= write("sonido.wav",99999,s)
