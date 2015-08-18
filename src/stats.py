import numpy
from scipy import stats
import matplotlib.pyplot as plt
from util import FileHandlers

class Correlation:
	def __init__(self):
		self.x_list = []
		self.y_list = []

	def _get_data(self, filename):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		txt_files = file_handlers.find_files(file_paths, 'txt')
		for txt_file in txt_files:
			if filename == file_handlers.get_file_name(txt_file):
				TXT = open(txt_file)
				data = TXT.readlines()
				TXT.close()
		return data

	def _build_dict(self, filename):
		values = {}
		data = self._get_data(filename)
		for line in data:
			resnum = line.split('\t')[0].strip()
			values[resnum] = line.split('\t')[1].strip()
		return values

	def _build_mapping(self, x_filename, y_filename):
		x_values = self._build_dict(x_filename)
		y_values = self._build_dict(y_filename)
		mapping = []
		for resnum in x_values:
			if resnum in y_values:
				#print resnum, x_values[resnum], y_values[resnum]
				mapping.append((x_values[resnum], y_values[resnum]))
		return mapping

	def linregress(self, x_filename, y_filename):
		mapping = self._build_mapping(x_filename, y_filename)
		x_values = []
		y_values = []
		for pair in mapping:
			if float(pair[0]) >= 0.35:
				x_values.append(float(pair[0]))
				y_values.append(float(pair[1]))
		slope, intercept, r_value, p_value, std_err = stats.linregress(x_values, y_values)
		print slope, intercept, r_value, p_value, std_err
		fig, ax = plt.subplots()
		fit = numpy.polyfit(x_values, y_values, deg=1)
		fit_fn = numpy.poly1d(fit)
		y = []
		for x in x_values:
			y.append(fit_fn(x))
		ax.plot(x_values, y, color='blue')
		ax.scatter(x_values, y_values)
		plt.show()














