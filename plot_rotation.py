import numpy as np
import matplotlib.pyplot as plt

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

CALIBRATION_FILE = "cal.csv"
DATA_FILE = 'tours.csv'
SEPARATOR = ';'


data = []



def integrate(t, y, origin = 0):
	out = [0]
	for i in range(1, max(len(t), len(y))):
		out.append(out[-1] + ((t[i] - t[i-1]) * (y[i] + y[i-1]) / 2))

	return out


def get_data(filename, offsets = None):
	data = []
	with open(filename) as file:
		line = file.readline()
		while line:
			data.append(line.split(SEPARATOR))
			line = file.readline()
	return [[float(r[i].replace(",", ".")) - (offsets[i] if offsets is not None else 0) for r in data[1:]] for i in range(len(data[0]))]

def detect_sat(d, std_dev, thres=20):
	count = 0
	c = [0]
	for i in range(1, len(d)):
		if abs(d[i - 1] - d[i]) < std_dev:
		#	print(count)
			count += 1
			if count > thres:
				return True
		else:
			count = 0
	return False


print("Getting calibration data")

cal_data = get_data(CALIBRATION_FILE)
offsets = [0]
std_dev = [0]

for data_set in cal_data[1:]:
	offsets.append(sum(data_set) / len(data_set))
	std_dev.append(np.std(data_set))

print("Offsets: ", offsets)
print("Standard deviation", std_dev)

print("Calibration done\n")

print("Detecting saturations")
data = get_data(DATA_FILE, offsets)

for i in range(1, len(data)):
	if detect_sat(data[i], std_dev[i]):
		print(bcolors.FAIL + "WARNING : Saturation detected" + bcolors.ENDC)

plt.show()

print("End of saturation detection\n")

print("Plotting")


for d in data[1:]:
	plt.plot(data[0], d)
	plt.plot(data[0], integrate(data[0], d))
plt.show()
