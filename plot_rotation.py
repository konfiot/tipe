import socket
from math import sqrt, atan2, asin, pi
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
DATA_FILE = 'rot90.csv'
SEPARATOR = ','


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
	with open(filename) as f:
		line = f.readline()
		while line:
			data.append(line.split(";"))
			line = f.readline()
	return [[float(r[i].replace(",", ".")) for r in data[1:]] for i in range(len(data[0]))]

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
	return sqrt(sum([x**2 for x in array]))

def normalize(array):
	return [x/norm(array) for x in array]

def get_raw_data(f):
	line, _ = f.recvfrom(8192)

	if not line:
		return False

	data = [float(x.replace(",", ".")) for x in line.split(SEPARATOR)]

	return [
		data[0], # Time
		data[6:9], # Rotation velocity
		data[2:5] # Linear acceleration
	]

print("Getting calibration data")

cal_data = get_cal_data(CALIBRATION_FILE)
offsets = []
std_dev = []

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


f = open(DATA_FILE, "r")
print("Data contains: " + f.readline())

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
s.bind(("", 5555))

plt.ion()
plt.show()
fig, ax = plt.subplots()
lines = ax.plot([],[],[],[],[],[])

ax.set_ylim(-pi, pi)

angles = [[], [], []]
t = []
c = 0
previous_data = [(TIME_INIT, QUAT_ROT_INIT, W_ARRAY_INIT)]
terminate = False

while not terminate:
	raw_data = get_raw_data(s)

	if not raw_data: # Escape the loop if there's no more data
		terminate = True
		continue

	for i in range(len(offsets)):
		raw_data[1][i] -= offsets[i]

	t.append(raw_data[0])


	p_t, p_quat, p_w = previous_data[-1] # Get the old data
	new_quat_rot = qmath.quaternion(norm(raw_data[1]) * (raw_data[0] - t[0] - p_t), normalize(raw_data[1])) # Calculate the new quaternion

	q = p_quat * new_quat_rot

	previous_data.append((raw_data[0] - t[0], q, raw_data[1])) # Add the new data to the old

	angles[0].append(atan2(2*q[0]*q[1] + q[2]*q[3], 1 - 2*(q[1]**2 + q[2]**2)))
	angles[1].append(asin(2*(q[0] * q[2] - q[3] * q[1])))
	angles[2].append(atan2(2*q[0]*q[3] + q[1]*q[2], 1 - 2*(q[2]**2 + q[3]**2)))

	print("%.5f" % angles[0][-1]) + "\t",
	print("%.5f" % angles[1][-1]) + "\t",
	print("%.5f" % angles[2][-1])

	c += 1

	if c%5 == 0:
		for i in range(len(lines)):
			lines[i].set_data(t, angles[i])
		plt.draw()
		ax.set_xlim(t[0], raw_data[0])

plt.ioff()
plt.show()


t = [v[0] for v in previous_data]
quat = [v[1] for v in previous_data]


