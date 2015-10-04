import math
import numpy as np
import matplotlib.pyplot as plt
import qmath

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


QUAT_ROT_INIT = qmath.quaternion(0, [0, 0, 0]) #Init the quaternion to the identity rotation (q = 1 + 0i + 0j + 0k = 1)
W_ARRAY_INIT = [0, 0, 0]
TIME_INIT = 0



def integrate(t, y, origin = 0):
	out = [0]
	for i in range(1, max(len(t), len(y))):
		out.append(out[-1] + ((t[i] - t[i-1]) * (y[i] + y[i-1]) / 2))

	return out


def get_cal_data(filename, offsets = None):
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
			count += 1
			if count > thres:
				return True
		else:
			count = 0
	return False

def norm(array):
	return math.sqrt(sum([x**2 for x in array]))

def normalize(array):
	return [x/norm(array) for x in array]

def get_raw_data(f):
	line = f.readline()

	if not line:
		return False


	data = [float(x.replace(",", ".")) for x in line.split(SEPARATOR)]

	return [
		data[0], # Time
		data[1:4], # Rotation velocity
		data[4:7] # Linear acceleration
	]

print("Getting calibration data")

cal_data = get_cal_data(CALIBRATION_FILE)
offsets = [0]
std_dev = [0]

for data_set in cal_data[1:]:
	offsets.append(sum(data_set) / len(data_set))
	std_dev.append(np.std(data_set))

print("Offsets: ", offsets)
print("Standard deviation", std_dev)

print("Calibration done\n")

#print("Detecting saturations")
#
#for i in range(1, len(data)):
#	if detect_sat(data[i], std_dev[i]):
#		print(bcolors.FAIL + "WARNING : Saturation detected" + bcolors.ENDC)
#
#print("End of saturation detection\n")

print("Integrating")

previous_data = [(TIME_INIT, QUAT_ROT_INIT, W_ARRAY_INIT)]
terminate = False

f = open(DATA_FILE, "r")
print("Data contains: " + f.readline())

while not terminate:
	raw_data = get_raw_data(f)

	if not raw_data: # Escape the loop if there's no more data
		terminate = True
		continue

	p_t, p_quat, p_w = previous_data[-1] # Get the old data
	new_quat_rot = qmath.quaternion(norm(raw_data[1]) * (raw_data[0] - p_t), normalize([raw_data[1][i] - offsets[i] for i in range(len(raw_data[1]))])) # Calculate the new quaternion

	previous_data.append((raw_data[0], p_quat * new_quat_rot, raw_data[1])) # Add the new data to the old



plt.show()
