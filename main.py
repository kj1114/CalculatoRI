#!/usr/bin/python
import numpy
import math
import csv
import unicodedata


from decimal import Decimal

# Function definition is here
def ErlangenWater( array ):
	"This prints a passed string into this function"
	# 	Wavelength should be given in microns

	#	Erlangen water
	T_star = 273.15;
	rho_star = 1000;
	lambda_star = 0.589;	 #micron

	a0 =0.244257733;
	a1 =9.74634476 * 10**-3;
	a2 = -3.73234996 * 10**-3;
	a3 =2.68678472*10**-4;
	a4 =1.58920570*10**-3;
	a5 =2.45934259*10**-3;
	a6 = 0.900704920;
	a7 =-1.66626219*10**-2;

	lambda_uv = 0.2292020;
	lambda_ir = 5.432937;

	lambda1 = numpy.array(array) / lambda_star

	T = 293/T_star;
	rho = 1;

	A = a0+a1*rho+a2*T+a3*lambda1**2*T+a4/lambda1**2+a5/(lambda1**2-lambda_uv**2)+a6/(lambda1**2-lambda_ir**2)+a7*rho**2
	n = numpy.array([math.sqrt(x) for x in ((2*A+1)/(1-A))]) - 0.001	# NB: substracted 0.001 to match experiments by Friebel and Olga

	return n.tolist();

def ErlangenWaterReference( array1, lambda_reference ):
   	"This function calculates ErlangenWater for a given reference value"

   	#This is where your reference lambda * 10^9 will be passed
   	#You need to fetch it, times by nm
   
   	# Find KK in Friebel data
   	# Find K in Sydoruk data
	if lambda_reference > 0.001:
		for index, entry in enumerate(array1):
			array1[index] = int(entry*10**9)

   	k = array1.index(lambda_reference)

   	for index, entry in enumerate(array1):
		array1[index] = entry*10**-9
   	# Calculate ErlangenWater at reference frequency
  	n_water_reference = [array1[k]]
   	n_water_reference = numpy.array(n_water_reference) * 10**6
   	n_water_reference = ErlangenWater(n_water_reference)


   	return k, n_water_reference

def SubstractiveKK(omega, Im_n, omega0, n0, k):
	"This function calculates SubstractiveKK"

	# omega0 and n0 are the reference values
	#k is the element number in omega and In_m corresponding to omega0

	#Remember that there are as many points as omega or Im_n_Hb

	#1. first value of frequency.
	#to calculate the principal value of the integral, exclude the first
	#value
	#secure a small interval around the integration point and around the
	#reference value.
	#To do so, connect two adjacent frequency points with a line, find the
	#line equation and step a little along this line
	
	ratio = 0.999 #setting the step

	N = len(Im_n)
	Re_n = []

	a = (Im_n[k+1]-Im_n[k])/(omega[k+1]-omega[k])
	b = (Im_n[k]*omega[k+1]-Im_n[k+1]*omega[k])/(omega[k+1]-omega[k])

	omega_plus = omega[k]*ratio
	Im_n_plus = a*omega_plus +b

	a = (Im_n[k]-Im_n[k-1])/(omega[k]-omega[k-1])
	b = (Im_n[k-1]*omega[k]-Im_n[k]*omega[k-1])/(omega[k]-omega[k-1])

	

	omega_minus = omega[k]/ratio
	Im_n_minus = a*omega_minus +b

	

	a = (Im_n[1]-Im_n[0])/(omega[1]-omega[0])
	b = (Im_n[0]*omega[1]-Im_n[1]*omega[0])/(omega[1]-omega[0])

	omega_close = omega[0]*ratio
	Im_n_close = a*omega_close +b;

	omega_dash = []
	Im_n_dash = []

	# the integration interval with secured circles
	indices = range(1, k)
	for x in indices:
		omega_dash.append(omega[x])
	omega_dash.append(omega_minus)
	omega_dash.insert(0, omega_close)\

	# the corresponding values of Im n
	for x in indices:
		Im_n_dash.append(Im_n[x])
	Im_n_dash.append(Im_n_minus)
	Im_n_dash.insert(0, Im_n_close)

	pod_integralom = []

	pod_integralom = (omega[0]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[0]**2)*((numpy.array(omega_dash))**2-omega0**2))
	I1 = numpy.trapz(pod_integralom, omega_dash)

	omega_dash = []
	Im_n_dash = []

	#second interval 
	indices = range(k+1, N)

	for x in indices:
		omega_dash.append(omega[x])
	omega_dash.insert(0, omega_plus)

	for x in indices:
		Im_n_dash.append(Im_n[x])
	Im_n_dash.insert(0, Im_n_plus)

	pod_integralom = (omega[0]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[0]**2)*((numpy.array(omega_dash))**2-omega0**2))
	I2 = numpy.trapz(pod_integralom, omega_dash)

	Re_temp = n0 - 2 / math.pi * (I1+I2)
	Re_temp = [Re_temp]
	
	for x in Re_temp:
		Re_n_first = x	# the integrals are negative because omega is descending

	#2. last value of frequency. Proceed as for the first one but choose
	#the intervals properly

	a = (Im_n[N-2]-Im_n[N-1])/(omega[N-2]-omega[N-1])
	b = (Im_n[N-1]*omega[N-2]-Im_n[N-2]*omega[N-1])/(omega[N-2]-omega[N-1])

	omega_close = omega[N-1]/ratio
	Im_n_close = a*omega_close +b

	omega_dash = []
	Im_n_dash = []

	indices = range(0, k)
	for x in indices:
		omega_dash.append(omega[x])
	omega_dash.append(omega_minus)

	for x in indices:
		Im_n_dash.append(Im_n[x])
	Im_n_dash.append(Im_n_minus)

	pod_integralom = (omega[N-1]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[N-1]**2)*((numpy.array(omega_dash))**2-omega0**2))
	I1 = numpy.trapz(pod_integralom, omega_dash)

	#second interval 

	omega_dash = []
	Im_n_dash = []

	indices = range(k+1, N-1)

	for x in indices:
		omega_dash.append(omega[x])
	omega_dash.append(omega_close)
	omega_dash.insert(0, omega_plus)

	for x in indices:
		Im_n_dash.append(Im_n[x])
	Im_n_dash.append(Im_n_close)
	Im_n_dash.insert(0, Im_n_plus)

	pod_integralom = (omega[N-1]**2-omega0**2) * numpy.array(omega_dash) * numpy.array(Im_n_dash)/((numpy.array(omega_dash)**2-omega[N-1]**2)*(numpy.array(omega_dash)**2-omega0**2))
	I2 = numpy.trapz(pod_integralom,omega_dash)

	Re_temp = n0 - 2 / math.pi * (I1+I2)
	Re_temp = [Re_temp]
	for x in Re_temp:
		Re_n_last = x

	#3. Now calculate for all the intermediate values of frequency
	#There are two loops, one if omega lies below the reference omega0,
	#and the other, if omega lies above omega0
	#The will be three intervals for integration in the two cases
	#omega[1]___omega-  ;  omega+___omega0- ;  omega0+___omega[N]
	#omega[1]___omega0-  ;  omega0+___omega- ;  omega+___omega[N]


	#Calculate for frequencies below omega0
	
	get_range = range(1,k)

	for i_ in get_range:

		#first interval

		a = (Im_n[i_-1]-Im_n[i_])/(omega[i_-1]-omega[i_])
		b = (Im_n[i_]*omega[i_-1]-Im_n[i_-1]*omega[i_])/(omega[i_-1]-omega[i_])


		omega_close = omega[i_]/ratio
		Im_n_close = a*omega_close +b

		indices = range(0, i_)
		omega_dash = []
		Im_n_dash = []

		for x in indices:
			omega_dash.append(omega[x])
		omega_dash.append(omega_close)

		for x in indices:
			Im_n_dash.append(Im_n[x])
		Im_n_dash.append(Im_n_close)

		pod_integralom = (omega[i_]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[i_]**2)*((numpy.array(omega_dash))**2-omega0**2))
		I1 = numpy.trapz(pod_integralom, omega_dash)

		#second interval 

		a = (Im_n[i_]-Im_n[i_+1])/(omega[i_]-omega[i_+1])
		b = (Im_n[i_+1]*omega[i_]-Im_n[i_]*omega[i_+1])/(omega[i_]-omega[i_+1])

		omega_close = omega[i_]*ratio
		Im_n_close = a*omega_close +b

		indices = range(i_+1, k)

		omega_dash = []
		Im_n_dash = []

		for x in indices:
			omega_dash.append(omega[x])
		omega_dash.append(omega_minus)
		omega_dash.insert(0, omega_close)

		for x in indices:
			Im_n_dash.append(Im_n[x])
		Im_n_dash.append(Im_n_minus)
		Im_n_dash.insert(0, Im_n_close)

		pod_integralom = (omega[i_]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[i_]**2)*((numpy.array(omega_dash))**2-omega0**2))
		I2 = numpy.trapz(pod_integralom, omega_dash)

		#third interval

		indices = range(k+1,N)
		omega_dash = []
		Im_n_dash = []

		for x in indices:
			omega_dash.append(omega[x])
		omega_dash.insert(0, omega_plus)

		for x in indices:
			Im_n_dash.append(Im_n[x])
		Im_n_dash.insert(0, Im_n_plus)

		pod_integralom = (omega[i_]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[i_]**2)*((numpy.array(omega_dash))**2-omega0**2))
		I3 = numpy.trapz(pod_integralom, omega_dash)

		Re_temp = n0 - 2 / math.pi * (I1+I2+I3)
		Re_temp = [Re_temp]
		for x in Re_temp:
			Re_n.append(x)

	for i_ in range (k+1, N-1):

		#First interval
		indices = range(0, k)

		omega_dash = []
		Im_n_dash = []

		for x in indices:
			omega_dash.append(omega[x])
		omega_dash.append(omega_minus)

		for x in indices:
			Im_n_dash.append(Im_n[x])
		Im_n_dash.append(Im_n_minus)

		pod_integralom = (omega[i_]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[i_]**2)*((numpy.array(omega_dash))**2-omega0**2))
		I1 = numpy.trapz(pod_integralom, omega_dash)

		#second interval 

		a = (Im_n[i_-1]-Im_n[i_])/(omega[i_-1]-omega[i_])
		b = (Im_n[i_]*omega[i_-1]-Im_n[i_-1]*omega[i_])/(omega[i_-1]-omega[i_])

		omega_close = omega[i_]/ratio
		Im_n_close = a*omega_close +b


		indices = range(k+1, i_)

		omega_dash = []
		Im_n_dash = []

		for x in indices:
			omega_dash.append(omega[x])
		omega_dash.append(omega_close)
		omega_dash.insert(0, omega_plus)

		for x in indices:
			Im_n_dash.append(Im_n[x])
		Im_n_dash.append(Im_n_close)
		Im_n_dash.insert(0, Im_n_plus)

		pod_integralom = (omega[i_]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[i_]**2)*((numpy.array(omega_dash))**2-omega0**2))
		I2 = numpy.trapz(pod_integralom, omega_dash)

		#third interval

		a = (Im_n[i_]-Im_n[i_+1])/(omega[i_]-omega[i_+1])
		b = (Im_n[i_+1]*omega[i_]-Im_n[i_]*omega[i_+1])/(omega[i_]-omega[i_+1])

		omega_close = omega[i_]*ratio
		Im_n_close = a*omega_close +b

		indices = range(i_+1, N)

		omega_dash = []
		Im_n_dash = []

		for x in indices:
			omega_dash.append(omega[x])
		omega_dash.insert(0, omega_close)

		for x in indices:
			Im_n_dash.append(Im_n[x])
		Im_n_dash.insert(0, Im_n_close)

		pod_integralom = (omega[i_]**2-omega0**2) * numpy.array(omega_dash)*numpy.array(Im_n_dash)/(((numpy.array(omega_dash))**2-omega[i_]**2)*((numpy.array(omega_dash))**2-omega0**2))
		I3 = numpy.trapz(pod_integralom, omega_dash)

		Re_temp = n0 - 2 / math.pi * (I1+I2+I3)
		Re_temp = [Re_temp]
		for x in Re_temp:
			Re_n.append(x)
	Re_n.insert(0, Re_n_first)
	Re_n.insert(k, n0) #at the reference frequency, there should be the reference value
	Re_n.append(Re_n_last)
	return numpy.array(Re_n)



def realPartRefractiveIndex(lambda_reference, check_type, concentration, Nreference, file1):

	
	n0_HbO2 = float(Nreference)				
	lambda_reference = int(lambda_reference)
	c = 299792458
	

	#request for is_checked missing
	#is_checked if checkbox ticked, otherwise false
	Re_n = []
	#Check if it is default setting 
	if file1 is None:
		concentration = float(concentration)
		#request for check_type missing
		#check_type if oxygeneted, otherwise deoxygeneted
		#For oxygeneted
		if (check_type!=0):
			wavelength = []
			Im_n_HbO2 = []
			e_HbO2 = []
			mu_a_HbO2 = []
			with open('prahl.txt') as f:
				for line in f:
					data = line.split()
					wavelength.append(float(data[0]))
					e_HbO2.append(float(data[1]))
				wavelength_cal = numpy.array(wavelength)*10**-9
				wavelength = wavelength_cal.tolist()

				omega = 2 * math.pi * c / wavelength_cal
				omega0 = 2 * math.pi * c / lambda_reference * 10**9

				mu_a_HbO2 = 2.303 * concentration / (64500) * 100  * numpy.array(e_HbO2)
				Im_n_HbO2 = mu_a_HbO2 * numpy.array(wavelength) /(4*math.pi)
				Im_n_HbO2 = Im_n_HbO2.tolist()

				[k, water_contribution] = ErlangenWaterReference(wavelength, lambda_reference)
				Re_n = SubstractiveKK(omega,Im_n_HbO2,omega0,n0_HbO2,k)+ErlangenWater(wavelength_cal*10**6)

		#For deoxygeneted
		else:
			wavelength = []
			Im_n_Hb = []
			e_Hb = []
			mu_a_Hb = []

			with open('prahl.txt') as f:
				for line in f:
					data = line.split()
					wavelength.append(float(data[0]))
					e_Hb.append(float(data[2]))
				wavelength_cal = numpy.array(wavelength)*10**-9
				wavelength = wavelength_cal.tolist()

				omega = 2 * math.pi * c / wavelength_cal
				omega0 = 2 * math.pi * c / lambda_reference * 10**9

				mu_a_Hb = 2.303*concentration/(64500)*100 * numpy.array(e_Hb)
				Im_n_Hb = mu_a_Hb * numpy.array(wavelength) /(4*math.pi)
				Im_n_Hb = Im_n_Hb.tolist()

				[k, water_contribution] = ErlangenWaterReference(wavelength, lambda_reference)	
				Re_n = SubstractiveKK(omega,Im_n_Hb,omega0,n0_HbO2,k)+ErlangenWater(wavelength_cal*10**6)

		result = [{'wavelength': wavelength_cal*10**9, 'Re_n': Re_n} for wavelength_cal, Re_n in zip(wavelength_cal, Re_n)]

		return result
	else:
		wavelength = []
		Im_n_HbO2 = []
		with file1 as f:
			for line in f:
				data = line.split()
				wavelength.append(float(data[0]))
				Im_n_HbO2.append(float(data[1]))
			if wavelength[0] > 0.001:
				wavelength_cal = numpy.array(wavelength)*10**-9
				wavelength = wavelength_cal.tolist()

				omega = 2 * math.pi * c / wavelength_cal
				omega0 = 2 * math.pi * c / lambda_reference * 10**9
				[k, water_contribution] = ErlangenWaterReference(wavelength, lambda_reference)				
				Re_n = SubstractiveKK(omega,Im_n_HbO2,omega0,n0_HbO2,k)+ErlangenWater(wavelength_cal*10**6)
			else:
				wavelength_cal = numpy.array(wavelength)

				omega = 2 * math.pi * c / wavelength_cal
				omega0 = 2 * math.pi * c / lambda_reference * 10**9
				[k, water_contribution] = ErlangenWaterReference(wavelength, lambda_reference)
				Re_n = SubstractiveKK(omega,Im_n_HbO2,omega0,n0_HbO2,k)+ErlangenWater(wavelength_cal*10**6)

		result = [{'wavelength': wavelength_cal*10**9, 'Re_n': Re_n} for wavelength_cal, Re_n in zip(wavelength_cal, Re_n)]
		return result



