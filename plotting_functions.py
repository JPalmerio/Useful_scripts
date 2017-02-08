from glob import glob
import os
try:
	from astropy.wcs import WCS
	from astropy.io import fits
except ImportError :
	print "Warning could not import astropy and pyfits. Cannot read fits files"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.transforms import blended_transform_factory
import numpy.random as rand
import scipy.stats as stats
import sys
import platform


""" Change this to you own home directory """
if platform.system() == 'Linux':
	homedir = '/nethome/palmerio/'
elif platform.system() == 'Darwin': 
	homedir = '/Users/palmerio/'
	

def where_value(value, array):
	"""
	Function that returns the indice of a value in a !!! sorted !!! array.
	Essentially, this is a dichotomy.

	Parameters :
	------------
	value : [float]
		Value to find in the array.

	array : [array]
		Array in which to look.

	Returns : 
	---------
	indice : [int]
		The indice where the value was found in the array.
	"""
	imin = 0
	imax = len(array)-1
	if(value <= array[imin] or value > array[imax]): 
		print "[Warning] in where_value : value is out of range of the array, indice returned : -1"
		indice = -1
	else :
		jj = 0 
		while True :
			j = int((imin + imax)/2)
			jj += 1
			if   (value >= array[j])  :
				imin = j+1
			elif (value < array[j-1]) :
				imax = j
			elif (jj > 1000):
				print "[Error] in where_value : looping over 1000 times, aborting."
				break
			else : 
				indice = j
				break
	return 	indice

def read_column(filename, column_nb, dtype=float, array=True, splitter=None, stripper=None):
	"""
	Function used to read ASCII (or fits) files.
	It will skip lines starting with '#', '!' or '%'.

	Parameters
	----------
	filename : [str]
		Name of the file containing the data

	column_nb: [int]
		Number of the column for the data (if fits file, 0 -> wavelength, 1 -> flux)

	dtype : [data-type]
		Type of the returned data. Default is float.

	array : [boolean]
		If True returns xdata as an array rather than a list. Default is True (arrays are faster)

	splitter : [str]
		String to use as a delimiter between columns. Default is None (uses default for str.split() which is a whitespace)

	stripper : [str]
		String to strip at the end of each line. Default is None (uses default for str.strip() which is a whitespace)

	Returns 
	-------
	xdata : [array/list]
		x data
	"""

	nan = False

	if filename[-5:] == '.fits' :
		xdata = read_fits3(filename)
		xdata = xdata[column_nb]
	else:
		xdata = []
		f = open(filename, 'r')
		for line in f :
			if len(line) != 0 :
				if line[0] != '#' and line[0] != '!' and line[0] != '%' :
					if stripper is not None:
						line = line.strip(stripper)
					else :
						line = line.strip()
					if splitter is not None:
						columns = line.split(splitter)
					else :
						columns = line.split()
					try :
						xdata.append(dtype(columns[column_nb]))
					except ValueError: 
						nan = True
						xdata.append(np.nan)
		if(array):
			xdata = np.asarray(xdata, dtype=dtype)
		f.close()
		if nan :
			print "[Warning] in read_column for %s. Could not convert string to %s, so added NaN." %(filename, dtype)

	return xdata

def read_data(filename, column_nb, err=True, dtype=float, splitter=None, stripper=None):
	"""
	Helper function that reads a table and returns the data, with errors and upper or lower limits.
	If err is True it assumes format of file is :     1) data     2) error plus     3) error minus 
	If there is no data, or if you ask for no errors, it fills with NaN 
	otherwise if the data exists but it doesn't find errors it'll set them to 0.
	Returns data as a list of numpy.ndarray with :
		data[0] = data (float)
		data[1] = plus error (float)
		data[2] = minus error (float)
		data[3] = upper limit (bool)
		data[4] = lower limit (bool)
	See help for read_column for more information on arguments.
	"""
	data = [[], [], [], [], []]
	f = open(filename, 'r')
	i = 1
	if err == True :
		for line in f:
			if line[0] != '#' and line[0] != '!' and line[0] != '%':
				if stripper is not None:
					line = line.strip(stripper)
				else :
					line = line.strip()
				if splitter is not None:
					columns = line.split(splitter)
				else :
					columns = line.split()
				try:
					# Check for upper limit
					if columns[column_nb][0] == '<':
						# Remove '<' from data
						columns[column_nb] = columns[column_nb][1:]
						# Append upper limit
						data[3].append(True)
						# Set errors to 0.0, and lower limit to false
						data[1].append(0.) # error plus
						data[2].append(0.) # error minus
						data[4].append(False)

					# Check for lower
					elif columns[column_nb][0] == '>':
						# Remove '>' from data
						columns[column_nb] = columns[column_nb][1:]
						# Append lower limit
						data[4].append(True)
						# Set errors to 0.0, and upper limit to false
						data[1].append(0.) # error plus
						data[2].append(0.) # error minus
						data[3].append(False)

					# if no lower or upper limits, append errors
					else : 
						data[3].append(False)
						data[4].append(False)
						try:
							data[1].append(float(columns[column_nb+1])) # error plus
						except IndexError:
							data[1].append(None) # error plus
						except ValueError:
							data[1].append(0.) # error plus
						try:
							data[2].append(np.abs(float(columns[column_nb+2]))) # error minus
						except IndexError:
							data[2].append(None) # error minus
						except ValueError:
							data[2].append(0.) # error minus
				
					# Append data
					data[0].append(dtype(columns[column_nb]))

				# If no data
				except IndexError:
					#print 'Warning : No data found for column %d, line %d in file %s. Input will be None.' %(column_nb, i, filename)
					data[0].append(None) # data
					data[1].append(None) # error plus
					data[2].append(None) # error minus
					data[3].append(None) # upper limit
					data[4].append(None) # lower limit
			i += 1
	else :
		for line in f:
			if line[0] != '#' and line[0] != '!' and line[0] != '%':
				if stripper is not None:
					line = line.strip(stripper)
				else :
					line = line.strip()
				if splitter is not None:
					columns = line.split(splitter)
				else :
					columns = line.split()
				try:
					# Check for upper limit
					if columns[column_nb][0] == '<':
						# Remove '<' from data
						columns[column_nb] = columns[column_nb][1:]
						# Append upper limit
						data[3].append(True)
						# Set lower limit to false
						data[4].append(False)

					# Check for lower
					elif columns[column_nb][0] == '>':
						# Remove '>' from data
						columns[column_nb] = columns[column_nb][1:]
						# Append lower limit
						data[4].append(True)
						# Set upper limit to false
						data[3].append(False)

					else :
						data[3].append(False)
						data[4].append(False)
					# Append data
					data[0].append(dtype(columns[column_nb]))
					data[1].append(None) # error plus
					data[2].append(None) # error minus


				except IndexError:
					#print 'Warning : No data found for column %d, line %d in file %s. Input will be None.' %(column_nb, i, filename)
					data[0].append(None) # data
					data[1].append(None) # error plus
					data[2].append(None) # error minus
					data[3].append(None) # upper limit
					data[4].append(None) # lower limit
	data[0] = np.array(data[0]).astype(dtype)
	data[1] = np.array(data[1]).astype(dtype)
	data[2] = np.array(data[2]).astype(dtype)
	data[3] = np.asarray(data[3])
	data[4] = np.asarray(data[4])
	return data

def read_overhead(filename, splitter=None, stripper=None):
	"""
		Helper function that reads the overhead that begins with '%' in an ASCII file.
		Returns a list of strings split with splitter.
	"""
	f = open(filename, 'r')
	overhead = []
	for line in f:
		if line[0]=='#' and line[1] == '%':
			line = line.strip('#')
			line = line.strip('%')
			if stripper is not None :
				line=line.strip(stripper)
			else :
				line=line.strip()
			if splitter is not None:
				line = line.split(splitter)
			else:
				line = line.split()
			for i in range(len(line)):
				if line[i] != '':
					overhead.append(line[i].strip())
	f.close()

	return overhead

def read_fits(filename):
    '''Read a UVES fits spectrum from the ESO pipeline

    Parameters
    ----------
    filename : [str]
    	Name of the fits file with the data

    Returns
    -------
    wavelength : [ndarray]
    	Wavelength (in Ang)

    flux : [ndarray]
    	Flux (in erg/s/cm**2)

    date_obs : [str]
    	Time of observation
    '''
    sp = fits.open(filename)
    header = sp[0].header

    wcs = WCS(header)
    #make index array
    index = np.arange(header['NAXIS1'])

    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    wavelength = wavelength.flatten()
    flux = sp[0].data

    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs

def read_fits3(name, ex=0):
    fitsfile1d=pyfits.open(name)
    data1d = fitsfile1d[ex].data
    n = len(data1d)
    lam=fitsfile1d[ex].header['CRVAL1']
    delta=fitsfile1d[ex].header['CDELT1']
    wavelength = np.zeros(n)
    for i in range(n):
        wavelength[i] = lam + i * delta
 
    return wavelength, data1d

def read_fits2(name, ex=1):
    fitsfile1d=pyfits.open(name)
    data1d = fitsfile1d[ex].data
    n = len(data1d)

    wavelength = [data1d[i][0] for i in range(n)]
    flux = [data1d[i][1] for i in range(n)]

    return wavelength, flux
  
def read_catalog(filename, ext=6):
	fitsfile = pyfits.open(filename)
	data1d = fitsfile[ext].data
	wavelength = np.zeros(len(data1d))
	flux = np.zeros(len(data1d))
	for i in range(len(data1d)):
		wavelength[i] = data1d[i][0]
		flux[i] = data1d[i][1]
	return wavelength, flux  

def smooth(data, sm_range):
	"""
	A function that averages the data over a number of points and returns the averaged data.
	Parameters
	----------
	data : [array]
		The data to be smoothed

	sm_range : [int]
		Range over which to average the data.

	Returns
	-------
	smoothed_data : [array]
		The data averaged over sm_range.

	"""
	n = len(data)%sm_range
	if (n == 0): n = len(data)/sm_range
	else : n = len(data)/sm_range + 1
	smoothed_data = np.zeros(n)
	if (sm_range > 1):		
		for i in range(n):
			if( (i+1)*sm_range <= len(data) ):
				smoothed_data[i] = np.mean(data[i*sm_range:(i+1)*sm_range])
			else :
				smoothed_data[i] = np.mean(data[i*sm_range:])
	else :
		smoothed_data = data
	return smoothed_data

def plot_spectrum_from_data(ax, wavelength, flux, sm_range=1, **kwargs):
	"""
	A helper function to plot a Spectrum from data

	Parameters
	----------
	ax : [Axes]
		The axes to draw to

	wavelength : [array]
		The x data (A)

	flux : [array]
		The y data (erg/cm2/s/A)

	**kwargs : [pointer]?
		Dictionary of kwargs to pass to ax.plot

	sm_range : [int]
		The number of pixels over which to smooth. Default is 1
	Returns
	-------
	spectrum : [list]
		List of artists added
	"""
	
	if (sm_range <= 1) :
		spectrum = ax.plot(wavelength, flux, **kwargs)
	else :
		smoothed_flux = smooth(flux, sm_range)
		smoothed_wvlg = smooth(wavelength, sm_range)
		spectrum = ax.plot(smoothed_wvlg, smoothed_flux, **kwargs)
	if wavelength[0] > 1000.:
		ax.set_xlabel(r'Wavelength ($\AA$)')
	else :
		ax.set_xlabel(r'Wavelength ($nm$)')

	ax.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ A$^{-1}$)')
	ax.grid(True)
	return spectrum

def plot_spectrum_from_file(ax, filename, column_nb_w=0, column_nb_f=1, sm_range=1, **kwargs):
	"""
	A helper function to plot a Spectrum from a file

	Parameters
	----------
	ax : [Axes]
		The axes to draw to

	filename : [str]
		Name of the file to plot

	column_nb_w : [int]
		Column number for the wavelength
	
	column_nb_f : [int]
		Column number for the flux

	**kwargs : [pointer]?
		Dictionary of kwargs to pass to ax.plot

	sm_range : [int]
			The number of pixels over which to smooth. Default is 1 (no smooth)

	Returns
	-------
	spectrum : [list]
		List of artists added
	"""
	wavelength = read_column(filename, column_nb_w)
	if np.min(wavelength) < 1000. :
		wavelength *= 10.
		print "[Info] In plot_spectrum_from_file, multiplying wavelength by 10 to convert to Angstroms."
	flux      = read_column(filename, column_nb_f)
	spectrum  = plot_spectrum_from_data(ax, wavelength, flux, sm_range=sm_range, **kwargs)
	return spectrum

def stack_fits(file_list, format='fits'):
	"""
	Check the assumption that wavelengths are the same for all files in the file_list.
	Output name is file_list[0] + 'stacked' in .fits format or .dat format, default is .fits

	"""

	n = len(file_list)
	#wavelengths = []
	fluxes = []
	for i in range(n):
		wavelength, flux = read_fits2(file_list[i])
		#wavelengths.append(wavelength)
		fluxes.append(flux)

	fluxes = np.asarray(fluxes)
	flux_stacked = np.mean(fluxes, axis=0)

	fitsfile1d = pyfits.open(file_list[0])
	head = fitsfile1d[0].header
	#	fitsfile1d[1].data=flux_stacked
	if format == 'fits':
		output_filename = file_list[0][:-5] + '_stacked.fits'
		pyfits.writeto(output_filename,flux_stacked, head)
	elif format == 'dat':
		output_filename = file_list[0][:-5] + '_stacked.dat'
		f = open(output_filename, 'w')
		for i in range(len(flux_stacked)):
			f.write(str(wavelength[i]) + '\t' + str(flux_stacked[i]) + '\n' )
		f.close()
	elif format == 'txt':
		output_filename = file_list[0][:-5] + '_stacked.txt'
		f = open(output_filename, 'w')
		for i in range(len(flux_stacked)):
			f.write(str(wavelength[i]) + '\t' + str(flux_stacked[i]) + '\n' )
		f.close()
	else :
		print "[Error] In stack_fits : wrong format. Please choose 'fits', 'dat', or 'txt'."
	return

def stack_fits2(file1,file2,file3, format=0):
	"""
	Check the assumption that wavelengths are the same for both files.
	Output name is file1 + 'stacked' in .fits format (0) or .dat format (1), default is .fits

	"""
	wavelength1, flux1 = read_fits2(file1)
	wavelength2, flux2 = read_fits2(file2)
	wavelength3, flux3 = read_fits2(file3)
	flux_stacked = (flux1 + flux2 + flux3)/3.0
	fitsfile1d = pyfits.open(file1)
	head = fitsfile1d[0].header
	#	fitsfile1d[1].data=flux_stacked
	if format == 0:
		output_filename = file1[:-5] + '_stacked.fits'
		pyfits.writeto(output_filename,flux_stacked, head)
	elif format == 1:
		output_filename = file1[:-5] + '_stacked.dat'
		f = open(output_filename, 'w')
		for i in range(len(flux_stacked)):
			f.write(str(wavelength1[i]) + '\t' + str(flux_stacked[i]) + '\n' )
		f.close()
	else :
		print "Error in stack_fits : wrong format"
	return

def spectral_gaussian(center, FWHM, flux, wavelength=None, continuum=0):
	"""	Function that returns a gaussian function.
		Note : FWHM = 2*sqrt(2*ln(2)) * sigma
		If wavelength is not specified, creates an arbitrary wavelength linspace.
		Otherwise it expects an array and will respect the sampling.
	"""
	sigma = FWHM / (2*np.sqrt(2*np.log(2)))
	height = flux / (np.sqrt(2*np.pi) * sigma)
	if wavelength == None:
		min_wvlg = center - 5*sigma
		max_wvlg = center + 5*sigma
		wavelength = np.linspace(min_wvlg, max_wvlg, 1000)
		print "Wavelength not specified, created linspace from %.1lf to %.1lf" %(min_wvlg, max_wvlg)
	gaussian = continuum + height * np.exp( -(wavelength - center)**2 / (2 * sigma**2) )

	return wavelength, gaussian

def convert_AB_to_Vega(AB_mag, band):
	"""
		Function that converts AB magnitudes to Vega magnitudes. 
		Comes from https://www.astro.umd.edu/~ssm/ASTR620/mags.html
		Accepted bands : V, B, Bj, R, I, g, r, i, u', g', r', i', z', Rc, Ic
	"""

	conversion = {'V':0.044, 'B':0.163, 'Bj':0.139, 'R':-0.055, 'I':-0.309, 'g':0.013,
				  'r':0.226, 'i':0.296, "u'":0.0, "g'":0.0, "r'":0.0, "i'":0.0, "z'":0.0,
				  'Rc':-0.117, 'Ic':-0.342}
	offset =  conversion.get(band, "Unvalid band")
	if offset == 'Unvalid band':
		print "Unvalid band in convert_AB_to_Vega.\nAccepted bands : V, B, Bj, R, I, g, r, i, u', g', r', i', z', Rc, Ic"
		Vega_mag = 0.0
	else:
		Vega_mag = AB_mag + offset

	return Vega_mag

def convert_Vega_to_AB(Vega_mag, band):
	"""
		Function that converts Vega magnitudes to Ab magnitudes. 
		Comes from https://www.astro.umd.edu/~ssm/ASTR620/mags.html
		Accepted bands : V, B, Bj, R, I, g, r, i, u', g', r', i', z', Rc, Ic
	"""

	conversion = {'V':0.044, 'B':0.163, 'Bj':0.139, 'R':-0.055, 'I':-0.309, 'g':0.013,
				  'r':0.226, 'i':0.296, "u'":0.0, "g'":0.0, "r'":0.0, "i'":0.0, "z'":0.0,
				  'Rc':-0.117, 'Ic':-0.342}
	offset =  conversion.get(band, "Unvalid band")
	if offset == 'Unvalid band':
		print "Unvalid band in convert_Vega_to_AB.\nAccepted bands : V, B, Bj, R, I, g, r, i, u', g', r', i', z', Rc, Ic"
		AB_mag = 0.0
	else:
		AB_mag = Vega_mag - offset

	return AB_mag

def convert_flux_to_AB_magnitude(flux, band, units='cgs'):
	""" 
		Function that converts a flux into AB magnitude.
		Assumes flux is in units of erg/s/cm2 (cgs) by default but can be changed to 'Jy'.
	"""
	if units == 'cgs':
		mag = -2.5 * np.log10(flux) - 48.6
	elif units == 'Jy':
		print "[Error] in convert_flux_to_AB_magnitude, Jy units not implemeted yet :D, sorry..."
		mag = 999.
	else :
		print "[Error] in convert_flux_to_AB_magnitude, invalid units. Please choose 'cgs', or 'Jy'. "
		mag = 999.
	return mag

def convert_AB_magnitude_to_flux(AB_mag, lam, units='cgs'):
	""" 
		Function that converts AB magnitude into flux.
		Lam is in angstroms.
		Returns flux in units of erg/s/cm2 (cgs) by default but can be changed to 'Jy'.
	"""

	if units == 'Jy' :
		flux = 3631*10.**(-0.4*AB_mag)
	elif units == 'cgs' :
		flux = 3631*10.**(-0.4*AB_mag) * 3 * 1e-5/lam**2
	else :
		print "[Error] in convert_AB_magnitude_to_flux, invalid units. Please choose 'cgs', or 'Jy'. "
	return flux

def plot_MCMC_path(filename, ax, log=[True, False], **kwargs):
	"""
	Helper function to easily plot the various MCMC paths.
	Filename can be a list of file names.
	Reads, x in column 0, y in column 1, and z in column 2.
	Takes the log of x if log[0] and y if log[1]. Default is [True, False]
	Replaces values >= 10**20 for z by np.nan.
	Returns a mappable for a colorbar.

	"""
	if type(filename) is list :
		for j in range(len(filename)):
			x     = read_data(filename[j], 0)
			if log[0] :
				x = np.log10(x)
			y     = read_data(filename[j], 1)
			if log[1] :
				y = np.log10(y)
			zchi2 = read_data(filename[j], 2) 
			for i in range(len(zchi2)):
				if zchi2[i] >= 10**20 :
					zchi2[i] = np.nan
			mask = np.isfinite(zchi2)
			cbn1 = ax.scatter(x[mask], y[mask], c=zchi2[mask], **kwargs)
			#ax.plot([x[0],x[-1]], [y[0],y[-1]], marker='D', markersize=50, color=['red', 'green'], ls=':', alpha=0.8)
			ax.plot([x[0],x[-1]], [y[0],y[-1]], marker=None, ls=':', color='k', zorder=190)
			ax.scatter([x[0],x[-1]], [y[0],y[-1]], marker='o', s=60, linewidths=0.5, c=['red', 'green'], zorder=200)


	else :
		x     = read_data(filename, 0)
		if log[0] :
			x = np.log10(x)
		y     = read_data(filename, 1)
		if log[1] :
			y = np.log10(y)
		zchi2 = read_data(filename, 2) 
		for i in range(len(zchi2)):
			if zchi2[i] >= 10**20 :
				zchi2[i] = np.nan
		mask = np.isfinite(zchi2)
		cbn1 = ax.scatter(x[mask], y[mask], c=zchi2[mask], **kwargs)

	return cbn1

def asym_gaussian_draw(mu, sigma1, sigma2, nb_draws=1000, precision=500, positive=False, ax=None, **kwargs):
	"""
	Function that draws randomly in a asymmetric gaussian distribution.
	in the form :   { exp( -(x-mu)**2/(2*sigma1**2) )     if x < mu
					{ exp( -(x-mu)**2/(2*sigma2**2) )     if x >= mu
	Also plots the distribution if ax is not None
	Returns an array of the drawings
	"""
	# Normalization
	#norm = 1. / ( np.sqrt(np.pi/2.) * (sigma1 + sigma2) )

	# errors need to be non zero so create artifical errors one billionth the value of mu (this can probably be better managed)
	if sigma1 == 0. :
		sigma1 = mu * 1e-9
	if sigma2 == 0. :
		sigma2 = mu * 1e-9

	# limits
	x_min = mu - 10.*sigma1	
	x_max = mu + 10.*sigma2	
	step = (x_max - x_min) / float(precision)

	# create Cumulative distribution 
	F_x = np.zeros(precision)
	p_x = np.zeros(precision)
	x = np.linspace(x_min, x_max, precision)
	# Initialization
	F_x[0] = 0.
	p_x[0] = gaussian(x[0], mu, sigma1)
	for i in range(precision-1):
		if x[i] < mu :
			p_x[i+1] = gaussian(x[i+1], mu, sigma1)
			p_x_a = gaussian(x[i], mu, sigma1)
			p_x_b = gaussian(x[i+1], mu, sigma1)
			F_x[i+1] = F_x[i] + step * (p_x_b + p_x_a)/2.
		elif x[i] >= mu :
			p_x[i+1] = gaussian(x[i+1], mu, sigma2)
			p_x_a = gaussian(x[i], mu, sigma2)
			p_x_b = gaussian(x[i+1], mu, sigma2)
			F_x[i+1] = F_x[i] + step * (p_x_b + p_x_a)/2.
	# normalize
	p_x /= F_x[-1]
	F_x /= F_x[-1] 

	# draw from generated distribution
	draw = np.zeros(nb_draws)
	for i in range(nb_draws):
		if positive:
			while draw[i] <= 0 :
				t = rand.random()
				j = where_value(t, F_x)
				draw[i] = x[j-1] + step*(t-F_x[j-1])/(F_x[j]-F_x[j-1])
		else:
			t = rand.random()
			j = where_value(t, F_x)
			draw[i] = x[j-1] + step*(t-F_x[j-1])/(F_x[j]-F_x[j-1])
	
	if ax is not None :
		ax.hist(draw, bins=30, normed=True, label=r'$\mu = %.2lf$' %mu + '\n' + r'$\sigma_- = %.2lf$' %sigma1 +'\n'+ r'$\sigma_+ = %.2lf$' %sigma2, **kwargs)
		ax.plot(x,p_x, color='k')
		ax.legend()

	return draw

def lin_to_log(x, x_errp, x_errm=None):
	"""
		Takes linear data with errors and converts to logscale with correct error propagation.
		If x_errm is not provided, errors are assumed symmetric.
		Returns : log_x, log_x_errp, log_x_errm
	"""
	if x_errm is None:
		x_errm = x_errp
	log_x_errp = np.log10(x + x_errp) - np.log10(x)
	log_x_errm = np.log10(x) - np.log10(x - x_errm)
	log_x      = np.log10(x)
	return log_x, log_x_errp, log_x_errm

def gaussian(x, mu, sigma):
	gaussian =  np.exp( -(x-mu)**2/(2*sigma**2) )
	return gaussian

def shifted_lines(z, units='A'):
	""" Prints the typical emission lines and their shifted wavelength in units or A or nm. """

	line_name = read_column(homedir+'Dropbox/Plotting_GUI/Lines/list_raie.dat', 0, dtype=str)
	line_wvlg = read_column(homedir+'Dropbox/Plotting_GUI/Lines/list_raie.dat', 1, dtype=float)

	for i in range(len(line_name)):
		if len(line_name[i]) < 9:
			line_name[i] += (9 - len(line_name[i]))*' '
		if units == 'A':
			print line_name[i] + '\t' + '%.3lf A' %(line_wvlg[i] * (1. + z)) + '\t'+' %.3lf'%line_wvlg[i]
		elif units == 'nm':
			print line_name[i] + '\t' + '%.3lf nm' %(line_wvlg[i] * (1. + z)/10.)+ '\t'+' %.3lf' %(line_wvlg[i]/10.)
		else :
			print " Wrong units, please choose 'A' or 'nm' "

def Lum_distance(z, units='Mpc'):
	"""
	Returns Luminosity Distance in Mpc by default but cm and Angstrom are also available.
	Uses Planck 2015 cosmology : Omega_M = 0.286, Omega_L = 0.714, H0 = 69.6 km/s/Mpc, no curvature.
	"""

	redshift_table = read_column(homedir+'Dropbox/python_codes/Observation/cosmo.dat', 0)
	D_L_table      = read_column(homedir+'Dropbox/python_codes/Observation/cosmo.dat', 1)
	i = where_value(z, redshift_table)

	D_L = D_L_table[i] 

	if units == 'cm':
		D_L *= 3.086*10**24
	elif units =='A':
		D_L *= 3.086*10**32

	return D_L

def star_formation_rate(flux, z, line='Ha', normalization='Chabrier'):
	"""
	Computes the SFR using H_alpha or [OII]3727 with different normalizations.
	Flux is interpreted as [erg/cm2/s].
	Note : if using OII, use combined flux of both lines.
	Also, if using for GRBs , need to add the relation from Savaglio and Kruhler.  
	Returns SFR in [M_sun/yr]. 
	"""

	norm = 0.

	if normalization == 'Chabrier':
		if line == 'Ha':
			norm = 4.6                # Chabrier
		elif line == 'OII':
			norm = 5.54 * 4.6/4.39    # Chabrier/Baldry 
		else:
			print "[Error] Wrong argument for 'line' in star_formation_rate. Use 'Ha' or 'OII'."
	elif normalization == 'Baldry':
		if line == 'Ha':
			norm = 4.39                # Baldry
		elif line == 'OII':
			norm = 5.54 * 4.39    
		else:
			print "[Error] Wrong argument for 'line' in star_formation_rate. Use 'Ha' or 'OII'."
	elif normalization == 'Kroupa':
		if line == 'Ha':
			norm = 5.37                # Kroupa
		elif line == 'OII':
			norm = 5.54 * 5.37/4.39    # Kroupa/Baldry
		else:
			print "[Error] Wrong argument for 'line' in star_formation_rate. Use 'Ha' or 'OII'."
	elif normalization == 'Salpeter':  # I haven't double checked this one 
		if line == 'Ha':
			norm = 7.9                # Salpeter
		elif line == 'OII':
			norm = 5.54 * 7.9/4.39    # Salpeter/Baldry
		else:
			print "[Error] Wrong argument for 'line' in star_formation_rate. Use 'Ha' or 'OII'."
	else:
		print "[Error] Wrong argument for 'normalization' in star_formation_rate.\nOptions are : ('Chabrier', 'Baldry', 'Kroupa', 'Salpeter')."

	D_L = Lum_distance(z, units='cm')
	Lum = 4. * np.pi * D_L**2 * flux 
	SFR = norm * 10**(-42) * Lum

	return SFR

def flux_to_lum(flux, z):
	"""
		Compute the Luminosity from a given flux in [erg/cm2/s] and redshift.
		Returns the luminosity in [erg/s].
	"""

	D_L = Lum_distance(z, 'cm')
	Lum = 4.0 * np.pi * D_L**2 * flux

	return Lum

def lum_to_flux(lum, z):
	"""
		Compute the flux from a given luminosity and redshift.
		Returns the flux in [erg/cm2/s].
	"""
	D_L = Lum_distance(z, 'cm')
	flux = lum / (4.0 * np.pi * D_L**2)

	return flux


def generate_colors(Nb_colors, seed=50):
	"""
		A helper function that returns a list of colors randomly generated of size Nb_colors.
		You can set the seed to avoid having changing colors everytime you replot.
	"""
	color = []
	np.random.seed(seed)
	for i in range(Nb_colors):
		color.append('#%06X' % np.random.randint(0, 0xFFFFFF))
	return color

class MC_var(object):
	"""
		Class to propagate errors using MC simulation.
		Needs a value, an error (wich is a list of two values : + and -), and a number of draws.
		You can also specify if you want you variable to be positive (False by default).
	"""

	def __init__(self, value, errorp, errorm=None, N=1000, pos=False):
		self.value = float(value)
		self.errorp = float(errorp)
		if errorm is None :
			self.errorm = self.errorp
		else :
			self.errorm = float(errorm)
		self.N = int(N)
		self.positive = pos
		#self.drawings = np.zeros(self.N)

	def __str__(self):
		return '%.3lf +/- (%.1lf, %.1lf)'%(self.value, self.errorp, self.errorm)

	def __add__(self, other):
		if type(other) == type(self):
			if other.N != self.N:
				newN = max(other.N, self.N)
			else :
				newN = self.N
	
			draw1 = self.asym_gaussian_draw(self.value, self.errorm, self.errorp, nb_draws=newN, positive=self.positive)
			draw2 = self.asym_gaussian_draw(other.value, other.errorm, other.errorp, nb_draws=newN, positive=other.positive)
			result = draw1 + draw2
			q = stats.mstats.mquantiles(result, prob=[0.16, 0.5, 0.84])
			value = q[1]
			errp  = q[2]-q[1]
			errm  = q[1]-q[0]
		else:
			value = other + self.value
			errp  = other + self.errorp
			errm  = other + self.errorm
			newN = self.N

		return MC_var(value, errp, errm, N=newN)

	def __sub__(self, other):
		if type(other) == type(self):
			if other.N != self.N:
				newN = max(other.N, self.N)
			else :
				newN = self.N
	
			draw1 = self.asym_gaussian_draw(self.value, self.errorm, self.errorp, nb_draws=newN, positive=self.positive)
			draw2 = self.asym_gaussian_draw(other.value, other.errorm, other.errorp, nb_draws=newN, positive=other.positive)
			result = draw1 - draw2
			q = stats.mstats.mquantiles(result, prob=[0.16, 0.5, 0.84])
			value = q[1]
			errp  = q[2]-q[1]
			errm  = q[1]-q[0]
		else:
			value = other - self.value
			errp  = other - self.errorp
			errm  = other - self.errorm
			newN = self.N

		return MC_var(value, errp, errm, N=newN)

	def __mul__(self, other):
		if type(other) == type(self):
			if other.N != self.N:
				newN = max(other.N, self.N)
			else :
				newN = self.N
	
			draw1 = self.asym_gaussian_draw(self.value, self.errorm, self.errorp, nb_draws=newN, positive=self.positive)
			draw2 = self.asym_gaussian_draw(other.value, other.errorm, other.errorp, nb_draws=newN, positive=other.positive)
			result = draw1 * draw2
			q = stats.mstats.mquantiles(result, prob=[0.16, 0.5, 0.84])
			value = q[1]
			errp  = q[2]-q[1]
			errm  = q[1]-q[0]

		else:
			value = other * self.value
			errp  = other * self.errorp
			errm  = other * self.errorm
			newN = self.N

		return MC_var(value, errp, errm, N=newN)



	def __div__(self, other):
		if type(other) == type(self):
			if other.N != self.N:
				newN = max(other.N, self.N)
			else :
				newN = self.N
	
			draw1 = self.asym_gaussian_draw(self.value, self.errorm, self.errorp, nb_draws=newN, positive=self.positive)
			draw2 = self.asym_gaussian_draw(other.value, other.errorm, other.errorp, nb_draws=newN, positive=other.positive)
			result = draw1 / draw2
			q = stats.mstats.mquantiles(result, prob=[0.16, 0.5, 0.84])
			value = q[1]
			errp  = q[2]-q[1]
			errm  = q[1]-q[0]

		else:
			value = self.value / other
			errp  = self.errorp / other
			errm  = self.errorm / other
			newN = self.N

		return MC_var(value, errp, errm, N=newN)

	def array(self):
		draw = self.asym_gaussian_draw(self.value, self.errorm, self.errorp, nb_draws=self.N, positive=self.positive)
		return draw

	def show(self, ax=None):
		if ax is None:
			ax = plt.subplot()

		self.asym_gaussian_draw(self.value, self.errorm, self.errorp, nb_draws=self.N, positive=self.positive, ax=ax)
		plt.draw()

	def asym_gaussian_draw(self, mu, sigma1, sigma2, nb_draws=5000, precision=1000, positive=False, ax=None, **kwargs):
		"""
		Function that draws randomly in a asymmetric gaussian distribution.
		in the form :   { exp( -(x-mu)**2/(2*sigma1**2) )     if x < mu
						{ exp( -(x-mu)**2/(2*sigma2**2) )     if x >= mu
		Also plots the distribution if ax is not None
		Returns an array of the drawings
		"""
		# Normalization
		#norm = 1. / ( np.sqrt(np.pi/2.) * (sigma1 + sigma2) )
	
		# errors need to be non zero so create artifical errors one billionth the value of mu (this can probably be better managed)
		if sigma1 == 0. :
			sigma1 = mu * 1e-9
		if sigma2 == 0. :
			sigma2 = mu * 1e-9
	
		# limits
		x_min = mu - 10.*sigma1	
		x_max = mu + 10.*sigma2	
		step = (x_max - x_min) / float(precision)
	
		# Cumulative distribution 
		F_x = np.zeros(precision)
		p_x = np.zeros(precision)
		x = np.linspace(x_min, x_max, precision)
		# Initialization
		F_x[0] = 0.
		p_x[0] = gaussian(x[0], mu, sigma1)
		for i in range(precision-1):
			if x[i] < mu :
				p_x[i+1] = gaussian(x[i+1], mu, sigma1)
				p_x_a = gaussian(x[i], mu, sigma1)
				p_x_b = gaussian(x[i+1], mu, sigma1)
				F_x[i+1] = F_x[i] + step * (p_x_b + p_x_a)/2.
			elif x[i] >= mu :
				p_x[i+1] = gaussian(x[i+1], mu, sigma2)
				p_x_a = gaussian(x[i], mu, sigma2)
				p_x_b = gaussian(x[i+1], mu, sigma2)
				F_x[i+1] = F_x[i] + step * (p_x_b + p_x_a)/2.
		# normalize
		p_x /= F_x[-1]
		F_x /= F_x[-1] 
	
		# draw from generated distribution
		draw = np.zeros(nb_draws)
		for i in range(nb_draws):
			if positive:
				while draw[i] <= 0 :
					t = rand.random()
					j = where_value(t, F_x)
					draw[i] = x[j-1] + step*(t-F_x[j-1])/(F_x[j]-F_x[j-1])
			else:
				t = rand.random()
				j = where_value(t, F_x)
				draw[i] = x[j-1] + step*(t-F_x[j-1])/(F_x[j]-F_x[j-1])
			if draw[i] == 0 :
				draw[i] = 1e-100
		
		if ax is not None :
			bins = np.arange(mu-5*sigma1, mu+5*sigma2, 15*step)
			ax.hist(draw, bins=bins, normed=True, label=r'$\mu = %.2lf$' %mu + '\n' + r'$\sigma_- = %.2lf$' %sigma1 +'\n'+ r'$\sigma_+ = %.2lf$' %sigma2, **kwargs)
			ax.plot(x,p_x, color='k')
			ax.legend(loc='best')
			ax.set_xlim([mu-5*sigma1, mu+5*sigma2])
	
		return draw
