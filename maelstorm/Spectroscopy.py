import numpy as n
import cataclysmic as cv

class Spectra(object):
	'''
	Class that contains all the data used in Doppler.
	JVHS - Dec 2013
	'''
	#import parameters as par
	#reload(par)

	def __init__(self, list):
		self.list = list
		self.files = {'name': [], 'phase': []}
		self.data = {'wave': [], 'flux': [], 'err': []}
		self.bin_data = {'wave': [], 'flux': [], 'err': []}
		self.Parameters = {'object': par.object,
							'bins': par.bins, 
							'porb': par.porb,
							'hjd0': par.hjd0,
							'gama': par.gama,
							'cenwave': par.lam0,
							'dellx' : par.dellx}
		self.author = 'Alymantara - University of Southampton - J.V.Hernandez@soton.ac.uk'
		self.medi = 3
		self.phase = n.arange(0+.5/self.Parameters['bins'],1,1./self.Parameters['bins']) 

	def Load_Spectra(self):
		'''
		Loads spectra and stores it in data dictionary
		Can do crazy stuff with the data

		'''
		self.data = {'wave': [], 'flux': [], 'err': []}
		self.files['name'],self.files['phase']=n.loadtxt(par.base_dir+'/'+self.list,dtype={'names': ('files', 'phase'),'formats': ('S12', 'f4')},unpack=True)
		for i in self.files['name']:
			print i
			w,f=n.loadtxt(par.base_dir+'/'+i,unpack=True)
			self.data['wave'].append(w),self.data['flux'].append(f)

	def Phase_Bin(self,plot_histo=True):
		'''
		Perfroms trail Spectra of the data loaded.

		You can modify the center wavelength and width of trail

		self.Parameters['cenwave']: Units of loaded spectra

		self.Parameters['dellx']: Pixels
		lc.Trail_Spectra(self.Parameters['cenwave']=6562.8,self.Parameters['dellx']=121)
		'''
		
		x,y=n.arange(0+.5/self.Parameters['bins'],2,1./self.Parameters['bins']),self.data['wave'][0]
		
		flag=0
		histo=[]

		for i in n.arange(self.Parameters['bins']):
		    tempfl=n.zeros(len(self.data['wave'][0]))
		    lolo=0
		    for j in n.arange(len(self.files['phase'])):
		        if i==0:
		            if self.files['phase'][j] < x[0] or self.files['phase'][j] >= 1-x[0]:
		                lolo=lolo+1
		                tempfl=cv.medfilt1(self.data['flux'][j],self.medi)+tempfl
		                #print 'on one ',1-bineq[0],' - ',bineq[0]
		        if self.files['phase'][j] >=x[i-1] and self.files['phase'][j]<x[i] and i!=0:
		            lolo=lolo+1
		            tempfl=cv.medfilt1(self.data['flux'][j],self.medi)+tempfl
		            #print 'after one>',bineq[i-1],' - ',bineq[i]
		    if lolo == 0:
		    	lolo =1
		    	medis=n.median(tempfl+1.0/lolo)
		    else:
		    	medis = n.median(tempfl/lolo)
		    #zvals[i]=tempfl/lolo/medis
		    self.bin_data['flux'].append(tempfl/lolo/medis)
		    #zvals[i+self.Parameters['bins']]=tempfl/lolo/medis
		    histo.append(lolo)
		

	def Trail_Spectra(self,cmaps=cm.Greys_r):
		import matplotlib.cm as cm
		import mynormalize
		import matplotlib.pyplot as plt
		ss=0
		for i in n.arange(len(self.data['wave'][0])-1):    
		    if self.Parameters['cenwave'] >= self.data['wave'][0][i]   and self.Parameters['cenwave'] <=self.data['wave'][0][i+1]:
		        ss=i
		x=self.data['wave'][0][ss-self.Parameters['dellx']:ss+self.Parameters['dellx']]
		zvals=[]
		for i in self.bin_data['flux']:
			zvals.append(i[ss-self.Parameters['dellx']:ss+self.Parameters['dellx']])
		zvals=n.array(zvals)
		fig = plt.figure(2,figsize=(8,10.5),facecolor='w')
		ax = fig.add_subplot(111)
		plt.clf()
		img = plt.imshow(self.bin_data['flux'],extent=(min(x), max(x),0, 2),interpolation='nearest', cmap=cmaps,aspect='auto')
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