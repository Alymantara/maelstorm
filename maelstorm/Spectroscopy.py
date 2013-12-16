import numpy as n
import cataclysmic as cv
import glob
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import resample as res
from pyfits import getval
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl

class Spectra(object):
	'''
	Class that contains all the data used in Doppler.
	JVHS - Dec 2013
	'''

	def __init__(self,param_file='parameters.py'):
		#self.list = list
		import re
		f=open(param_file)
		lines=f.readlines()
		f.close()

		self.Parameters = {'object': str(re.search('%s(.*)%s' % ('"', '"'), filter(lambda x: 'object' in x,lines)[0]).group(1)),
			'bins': int(re.search('%s(.*)%s' % ('=', '#'), filter(lambda x: 'bins' in x,lines)[0].replace('\t',' ')).group(1)), 
			'porb': float(re.search('%s(.*)%s' % ('=', '#'), filter(lambda x: 'porb' in x,lines)[0].replace('\t',' ')).group(1)),
			'hjd0': float(re.search('%s(.*)%s' % ('=', '#'), filter(lambda x: 'hjd0' in x,lines)[0].replace('\t',' ')).group(1)),
			'gama': float(re.search('%s(.*)%s' % ('=', '#' ), filter(lambda x: 'gama' in x,lines)[0].replace('\t',' ')).group(1)),
			'cenwave': float(re.search('%s(.*)%s' % ('=', '#'), filter(lambda x: 'lam0' in x,lines)[0].replace('\t',' ')).group(1)),
			'dellx' : int(re.search('%s(.*)%s' % ('=', '#'), filter(lambda x: 'dellx' in x,lines)[0].replace('\t',' ')).group(1)),
			'base_dir':str(re.search('%s(.*)%s' % ('"', '"'), filter(lambda x: 'base_dir' in x,lines)[0].replace('\t',' ')).group(1))}
		self.phase = n.arange(0+.5/self.Parameters['bins'],1,1./self.Parameters['bins']) 

		self.files = {'name': [], 'phase': []}
		self.data = {'wave': [], 'flux': [], 'err': []}
		self.bin_data = {'wave': [], 'flux': [], 'err': []}
		self.author = 'Alymantara - University of Southampton - J.V.Hernandez@soton.ac.uk'
		self.medi = 3

		 

	def Load_Spectra(self):
		'''
		Loads spectra and stores it in data dictionary
		Can do crazy stuff with the data

		'''
		self.data = {'wave': [], 'flux': [], 'err': []}
		self.files['name'],self.files['phase']=n.loadtxt(Parameters.base_dir+'/'+self.list,dtype={'names': ('files', 'phase'),'formats': ('S12', 'f4')},unpack=True)
		for i in self.files['name']:
			print i
			w,f=n.loadtxt(Parameters.base_dir+'/'+i,unpack=True)
			self.data['wave'].append(w),self.data['flux'].append(f)


	def Read_spectra(self,type='xshooter',output='input_nir',median_filter=False,resample=False,min_f=0,max_f=24):
		'''
		Reads spectra directly from fits file. Calculates Phase for every image.

		-----------

		output:  name for text file and *.npy file for easy retrieval after reading

		median_filter :  Applies a median filter to each spectra of N pixels.

		resample: 			Rebin all spectra to first linear dispersion encountered.

		'''
		self.data = {'wave': [], 'flux': [], 'err': [],'phase': []}

		if type == 'text':
			self.files['name'],self.files['phase']=n.loadtxt(Parameters.base_dir+'/'+self.list,dtype={'names': ('files', 'phase'),'formats': ('S12', 'f4')},unpack=True)
			for i in self.files['name']:
				print 'Loading: ', i.split('/')[-1][:-5]
				w,f=n.loadtxt(Parameters.base_dir+'/'+i,unpack=True)
				self.data['wave'].append(w),self.data['flux'].append(f)


		if type == 'xshooter':
			arm='NIR'
			files=glob.glob('/Users/juan/astro/SDSS1433/spectroscopy/'+arm+'/*'+arm+'.fits')

			
			wave,flux,name,mjd,ra,decl,phase,delphase,files1,hjd=[],[],[],[],[],[],[],[],[],[]
			wt,ft=cv.read_xshooter(files[0])
			f=open(output+'.txt','w')
			for i in files[min_f:max_f]:
			    wav,flu=cv.read_xshooter(i)

			    if resample:
			        flu=res.rebin2(wav,flu,wt)
			        flu[n.isnan(flu)] = 0.0
			    t1,t2,t3,t4,t5=getval(i,'OBJECT'),getval(i,'MJD-OBS'),getval(i,'RA'),getval(i,'DEC'),getval(i,'EXPTIME')
			    if median_filter:
			    	wave.append(wav),flux.append(cv.medfilt1(flu,17)),name.append(t1),mjd.append(t2),ra.append(t3),decl.append(t4),delphase.append(t5/3600/24/self.Parameters['porb'])
			    else:
			    	wave.append(wav),flux.append(flu),name.append(t1),mjd.append(t2),ra.append(t3),decl.append(t4),delphase.append(t5/3600/24/self.Parameters['porb'])

			    hjd.append(pyasl.helio_jd(t2+0.5,218.324772,10.19017))
			    phase.append(cv.phase(hjd[-1],self.Parameters['hjd0'],self.Parameters['porb']))
			    files1.append(i)
			    print 'Loading: '+i.split('/')[-1][:-5],' ' , t2,' ' ,hjd[-1],' ' ,phase[-1]
			    print >>f, i.split('/')[-1][:-5]+'.txt' ,phase[-1]
			    self.data['wave'].append(wav*10.0),self.data['flux'].append(flu),self.data['phase'].append(phase[-1])
			f.close()
			n.save(output,[wave,flux,name,hjd,phase,delphase,files1])

	def Phase_Bin(self,plot_histo=True):
		'''
		Perfroms trail Spectra of the data loaded.

		You can modify the center wavelength and width of trail

		self.Parameters['cenwave']: Units of loaded spectra

		self.Parameters['dellx']: Pixels
		lc.Trail_Spectra(self.Parameters['cenwave']=6562.8,self.Parameters['dellx']=121)
		'''
		self.bin_data['flux']=[]
		#x,y=n.arange(0+.5/self.Parameters['bins'],2,1./self.Parameters['bins']),self.data['wave'][0]
		bineq=n.arange(0+.5/self.Parameters['bins'],2,1./self.Parameters['bins'])

		flag=0
		histo=[]

		for i in n.arange(self.Parameters['bins']):
		    tempfl=n.zeros(len(self.data['wave'][0]))
		    lolo=0
		    cv.Printer( 'Binning '+str(float(i)/self.Parameters['bins']))
		    for j in n.arange(len(self.data['phase'])):
		        if i==0:
		            if self.data['phase'][j] < bineq[0] or self.data['phase'][j] >= 1-bineq[0]:
		                lolo=lolo+1
		                tempfl=cv.medfilt1(self.data['flux'][j],self.medi)+tempfl
		                #print 'on one ',1-bineq[0],' - ',bineq[0]
		        if self.data['phase'][j] >=bineq[i-1] and self.data['phase'][j]<bineq[i] and i!=0:
		            lolo=lolo+1
		            tempfl=cv.medfilt1(self.data['flux'][j],self.medi)+tempfl
		            #print 'after one>',bineq[i-1],' - ',bineq[i]
		    if lolo == 0:
		    	lolo =1
		    	medis=n.median((tempfl+1.0)/lolo)
		    else:
		    	medis = n.median(tempfl/lolo)
		    #zvals[i]=tempfl/lolo/medis
		    self.bin_data['flux'].append(tempfl)
		    #zvals[i+self.Parameters['bins']]=tempfl/lolo/medis
		    histo.append(lolo)
		if plot_histo:
			fig = plt.figure(figsize=(3,4))
			plt.step(bineq[:self.Parameters['bins']],histo,'b',where='mid')
			plt.step([i+1.0 for i in bineq[:self.Parameters['bins']]],histo,'b',where='pre')
			plt.step([i-1.0 for i in bineq[:self.Parameters['bins']]],histo,'b',where='mid')
			plt.axis([0,1.0,0,max(histo)+1])

	def Trail_Spectra(self,cmaps=cm.Greys_r):
		import matplotlib.cm as cm
		import mynormalize
		import matplotlib.pyplot as plt
		ss=0
		for i in n.arange(len(self.data['wave'][0])-1):    
		    if self.Parameters['cenwave'] >= self.data['wave'][0][i]   and self.Parameters['cenwave'] <=self.data['wave'][0][i+1]:
		        ss=i
		#x,y=n.arange(0+.5/self.Parameters['bins'],2,1./self.Parameters['bins']),self.data['wave'][0]
		x=self.data['wave'][0][ss-self.Parameters['dellx']:ss+self.Parameters['dellx']]
		zvals=[]
		for i in self.bin_data['flux']:
			zvals.append(i[ss-self.Parameters['dellx']:ss+self.Parameters['dellx']]/median(i[ss-self.Parameters['dellx']:ss+self.Parameters['dellx']]))
		zvals=n.array(zvals)
		fig = plt.figure(2,figsize=(8,10.5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.clf()
		img = plt.imshow(zvals,extent=(min(x), max(x),0, 2),interpolation='nearest', cmap=cmaps,aspect='auto')
		cbar = plt.colorbar(format='%05.2f')
		cbar.set_label('Arbitrary Flux')
		cbar.set_norm(mynormalize.MyNormalize(vmin=zvals.min(),vmax=zvals.max(),stretch='linear'))
		cbar = DraggableColorbar(cbar,img)
		cbar.connect()
		'''
		if plot_histo == True:
			fig = plt.figure(11,figsize=(2,3))
			plt.clf()
			plt.step(bineq[:bins],histo,'b',where='mid')
			plt.step([i+1.0 for i in bineq[:bins]],histo,'b',where='pre')
			plt.step([i-1.0 for i in bineq[:bins]],histo,'b',where='mid')
			plt.axis([0,1.0,0,max(histo)+1])
			plt.title('Spectra per bin. Total bins: '+str(bins))
			plt.ylabel('Number of spectra per bin')
			plt.xlabel('Phase Bin')

		# Plots the Trailed spectra
		fig = plt.figure(2,figsize=(8,10.5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.clf()
		if xaxis is 'wave':
		    img = plt.imshow(zvals,extent=(y.min(), y.max(),x.min()-0.5/bins, x.max()+0.5/bins),interpolation='nearest', cmap=cmaps,aspect='auto',origin='lower',vmin=minim,vmax=maxim)
		    plt.xlabel('Wavelength, $nm$')
		    limits=[y.min(), y.max(),x.min(), x.max()]
		    plt.axvline(x=lam,linestyle='--',color='white')
		if xaxis is 'vel':
		    velo=cv.redshift1(saved[0][0][ss-dell:ss+dell],lam)
		    img = plt.imshow(zvals,interpolation='nearest', cmap=cmaps,aspect='auto',origin='lower',extent=(min(velo), max(velo),x.min()-0.5/bins, x.max()+0.5/bins),vmin=minim,vmax=maxim)
		    plt.xlabel('Velocity, [km s$^{-1}$]')
		    limits=[min(velo), max(velo),x.min()-0.5/bins, x.max()+0.5/bins]
		    plt.axvline(x=0.0,linestyle='--',color='white')
		#extent=(min(velo), max(velo),x.min(), x.max())  ,vmin=0.2,vmax=1.3
		#cbar = fig.colorbar(img, ticks=[minim, minim +(maxim-minim)/5.,minim +(maxim-minim)*2/5.,minim +(maxim-minim)*3/5.,minim +(maxim-minim)*4/5.,maxim ])
		cbar = plt.colorbar(format='%05.2f')
		cbar.set_label('Arbitrary Flux')
		cbar.set_norm(mynormalize.MyNormalize(vmin=zvals.min(),vmax=zvals.max(),stretch='linear'))
		cbar = DraggableColorbar(cbar,img)
		cbar.connect()

		plt.show()
		plt.ylabel('Orbital Phase')
		plt.suptitle('CataPy Trailed Spectra')
		plt.title('Object: '+saved[2][0]+'\n Line: '+line_lbl+' - $\lambda_0=$'+str(lam)+' nm, $\Delta$pix='+str(dell)+', $\Delta\phi=$'+str(cv.trunc(1./bins,3))+'\n')
		phis=bineq-1./bins
		'''

	'''
    def foldspec(rebin=True,nbins=None):
   
    
    	nbins=None
    	inputs=n.loadtxt(par.base_dir+'/'+par.list,dtype={'names': ('files', 'phase'),'formats': ('S12', 'f4')})

    # Check 1st spectrum and get wavelength to interpolate
    	w1st=n.loadtxt(par.base_dir+'/'+inputs[0]['files'],unpack=True)
    	if nbins==None:
        	nbins=len(inputs)*1.5    #By default 
    	wave,flux=[],[]
    	for i in inputs:
        	print i['files']+'  '+str(i['phase'])
        	w,f=n.loadtxt(par.base_dir+'/'+i['files'],unpack=True)
        	#Missing interpolation in the case of different dispersion in wavelength
        	wave.append(w),flux.append(f)
    	delp=1.0/nbins
    	pha=n.arange(0,1,delp)    


    	return(wave,flux,pha,inputs['phase'],inputs['phase'],delp,inputs['files'])
    '''


class DraggableColorbar(object):
	def __init__(self, cbar, mappable):
		self.cbar = cbar
		self.mappable = mappable
		self.press = None
		self.cycle = sorted([i for i in dir(plt.cm) if hasattr(getattr(plt.cm,i),'N')])
		self.index = self.cycle.index(cbar.get_cmap().name)

	def connect(self):
		"""connect to all the events we need"""
		self.cidpress = self.cbar.patch.figure.canvas.mpl_connect(
			'button_press_event', self.on_press)
		self.cidrelease = self.cbar.patch.figure.canvas.mpl_connect(
			'button_release_event', self.on_release)
		self.cidmotion = self.cbar.patch.figure.canvas.mpl_connect(
			'motion_notify_event', self.on_motion)
		self.keypress = self.cbar.patch.figure.canvas.mpl_connect(
			'key_press_event', self.key_press)

	def on_press(self, event):
		"""on button press we will see if the mouse is over us and store some data"""
		if event.inaxes != self.cbar.ax: return
		self.press = event.x, event.y

	def key_press(self, event):
		if event.key=='down':
			self.index += 1
		elif event.key=='up':
			self.index -= 1
		if self.index<0:
			self.index = len(self.cycle)
		elif self.index>=len(self.cycle):
			self.index = 0
		cmap = self.cycle[self.index]
		self.cbar.set_cmap(cmap)
		self.cbar.draw_all()
		self.mappable.set_cmap(cmap)
		#self.mappable.get_axes().set_title(cmap)
		self.cbar.patch.figure.canvas.draw()

	def on_motion(self, event):
		'on motion we will move the rect if the mouse is over us'
		if self.press is None: return
		if event.inaxes != self.cbar.ax: return
		xprev, yprev = self.press
		dx = event.x - xprev
		dy = event.y - yprev
		self.press = event.x,event.y
		#print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
		scale = self.cbar.norm.vmax - self.cbar.norm.vmin
		perc = 0.03
		if event.button==1:
			self.cbar.norm.vmin -= (perc*scale)*n.sign(dy)
			self.cbar.norm.vmax -= (perc*scale)*n.sign(dy)
		elif event.button==3:
			self.cbar.norm.vmin -= (perc*scale)*n.sign(dy)
			self.cbar.norm.vmax += (perc*scale)*n.sign(dy)
		self.cbar.draw_all()
		self.mappable.set_norm(self.cbar.norm)
		self.cbar.patch.figure.canvas.draw()


	def on_release(self, event):
		"""on release we reset the press data"""
		self.press = None
		self.mappable.set_norm(self.cbar.norm)
		self.cbar.patch.figure.canvas.draw()

	def disconnect(self):
		"""disconnect all the stored connection ids"""
		self.cbar.patch.figure.canvas.mpl_disconnect(self.cidpress)
		self.cbar.patch.figure.canvas.mpl_disconnect(self.cidrelease)
		self.cbar.patch.figure.canvas.mpl_disconnect(self.cidmotion)