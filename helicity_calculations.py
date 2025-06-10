from pycwt import wavelet
import numpy as np

def helicity(Bbx,Bby,Bbz,Vx,Vy,Vz,Np,dt,dj,magnetosheath=False):

	# ALFVEN SPEED AND ELSASSER VARIABLES
	B_vec = np.stack((Bbx, Bby, Bbz))
	V_vec = np.stack((Vx, Vy, Vz))
	b = B_vec*(1e-9)/np.sqrt(4*np.pi*Np*(1.67e-27)*(1e-1))*(1e-3) #km/s
	Zplus = V_vec+b #km/s
	Zminus = V_vec-b #km/s

	# WAVELET TRANSFORMS
	#dt = 60 #seconds
	#dj = 0.125
	s0 = 2*dt
	mother = wavelet.Morlet(6) #2.5 Zhao

	Wx_kt, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Bbx, dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	Wy_kt, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Bby, dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	Wz_kt, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Bbz, dt, dj=dj, J=-1, s0=s0, wavelet=mother)

	wave_Vx, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Vx, dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Vy, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Vy, dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Vz, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Vz, dt, dj=dj, J=-1, s0=s0, wavelet=mother)

	wave_Zxplus, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Zplus[0], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Zyplus, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Zplus[1], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Zzplus, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Zplus[2], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Zxminus, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Zminus[0], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Zyminus, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Zminus[1], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	wave_Zzminus, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(Zminus[2], dt, dj=dj, J=-1, s0=s0, wavelet=mother)

	Wx_kt2, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(b[0], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	Wy_kt2, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(b[1], dt, dj=dj, J=-1, s0=s0, wavelet=mother)
	Wz_kt2, scale, freqs, coi, fft, fftfreqs = wavelet.cwt(b[2], dt, dj=dj, J=-1, s0=s0, wavelet=mother)

	# ELSASSER TRANSFORMS
	Wbp_kt = np.abs(wave_Zxplus)**2 + np.abs(wave_Zyplus)**2 + np.abs(wave_Zzplus)**2
	Wbm_kt = np.abs(wave_Zxminus)**2 + np.abs(wave_Zyminus)**2 + np.abs(wave_Zzminus)**2

	# KINETIC AND MAGNETIC ENERGY
	Wk_kt = np.abs(wave_Vx)**2 + np.abs(wave_Vy)**2 + np.abs(wave_Vz)**2
	Wm_kt = np.abs(Wx_kt2)**2 + np.abs(Wy_kt2)**2 + np.abs(Wz_kt2)**2

	# MAGNETIC HELICITY
	Hm = 2*(Wy_kt*np.conj(Wz_kt)).imag	# same as 2(-np.conj(Wy_kt)*Wz_kt).imag
	sig_m = Hm/(np.abs(Wx_kt)**2 + np.abs(Wy_kt)**2 + np.abs(Wz_kt)**2)

	if magnetosheath:
		N = len(Vx)
		localVx = np.empty((scale.size,N))
		localVy = np.empty((scale.size,N))
		localVz = np.empty((scale.size,N))
		for s in range(scale.size):
			for t in range(N):
				localVx[s,t] = np.sum(Vx*np.exp(-((t*dt - np.arange(0,N*dt,step=dt))**2)/(2*(scale[s]**2))))
				localVy[s,t] = np.sum(Vy*np.exp(-((t*dt - np.arange(0,N*dt,step=dt))**2)/(2*(scale[s]**2))))
				localVz[s,t] = np.sum(Vz*np.exp(-((t*dt - np.arange(0,N*dt,step=dt))**2)/(2*(scale[s]**2))))

		Hm = 2*(localVx*(-np.conj(Wy_kt)*Wz_kt).imag ) #+ localVy*(np.conj(Wz_kt)*Wx_kt).imag + localVz*(np.conj(Wx_kt)*Wy_kt).imag)
		sig_m = Hm/np.sqrt(localVx**2 + localVy**2 + localVz**2)/(np.abs(Wx_kt)**2 + np.abs(Wy_kt)**2 + np.abs(Wz_kt)**2)

	# CROSS HELICITY
	Hc = Wbp_kt-Wbm_kt
	sig_c = Hc/(Wbp_kt+Wbm_kt)

	# RESIDUAL ENERGY
	Er = Wk_kt-Wm_kt
	sig_r = Er/(Wk_kt+Wm_kt)

	return [scale, coi, sig_m, sig_c, sig_r]
