import numpy as np
import matplotlib.pyplot as plt
import csv

#PATH1 = '/Users/alihanks/Google Drive/NQUAKE_analysis/PERM/PERM_data/lbnl_sensor_60.csv'
PATH1 = '/Users/alihanks/k40_test_2019-02-06_D3S.csv'

def make_int(lst):
	'''
	Makes all entries of a list an integer
	'''
	y = []
	for i in lst:
		y.append(int(i))
	return y

def make_array(lst):
	'''
	Makes list into an array. Also splices out the irrelevant stuff
	for a spectra
	'''
	y = np.asarray(make_int(lst[12:]))
	return y

def main_potassium(number, n=1, lower_limit=270, upper_limit=292):
	'''
	Main Function.
	Number is the number of spectras to go through (Right now the total number is being used)
	n is the number of hours that each spectra is integrated over. Default is 1.
	It was experimentally determined that the channel for the K-40 centroid lies between 270 and 292 for the current
	peak finder. If a centroid falls outside that range (the lower_limit and upper_limit), then it is put in the anomaly list
	and it is plotted.
	In order to plot individual spectra, tune the lower_limit and upper_limit (these only plot spectra outside the range)
	'''

	entries = 300*n
	anomaly = []
	days = (24/n)
	i = 0
	counter = 0
	indexes = []
	day = 1
	while i < number:
		if counter < days:
			first_integration = rows[(i*entries)+1:((i+1)*entries)+1]
			array_lst = []
			for j in first_integration:
				array_lst.append(make_array(j))

			integrated = sum(array_lst)

			fig, ax = plt.subplots()
			fig.patch.set_facecolor('white')
			plt.title('PERM Spectra')
			plt.xlabel('channels')
			plt.ylabel('counts')
			plt.xlim(1,1000)
			#plt.ylim()
			x = range(0,len(integrated))
			ax.plot(x, integrated, 'bo-', label="CPM")
			ax.errorbar(x, integrated, yerr=np.sqrt(integrated), fmt='bo', ecolor='b')
			plt.yscale('log')
			plt.show()
			i+=1
			counter +=1
		else:
			#plt.title('1460 Centroid versus Time for Day {}'.format(day))
			#plt.xlabel('hours')
			#plt.ylabel('1460 Centroid')
			#plt.plot(indexes, 'ro')
			#plt.ylim(260, 300)
			#plt.show()
			print('plotted', day)
			counter = 0
			indexes = []
			day += 1
	#plt.title('1460 Centroid versus Time for Day {}'.format(day))
	#plt.xlabel('Hour of the Day')
	#plt.ylabel('1460 Centroid')
	#plt.plot(indexes, 'ro')
	#plt.ylim(260, 300)
	#plt.show()
	print('plotted', day)
	counter = 0
	indexes = []
	day += 1
	if anomaly:
		print(anomaly)
	else:
		print('There are no anomalies')


if __name__ == '__main__':
	with open(PATH1) as f:
	    reader = csv.reader(f)
	    rows = [r for r in reader]

	print('This data is taken from the {} csv'.format(PATH1))
	main_potassium(len(rows), n=1, lower_limit=0, upper_limit=4096)
