#!/usr/bin/python3

import socket
from math import sqrt, atan2, asin, pi
import numpy as np
import matplotlib.pyplot as plt
from mathutils import Quaternion
import math

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

CALIBRATION_FILE_GYROS = "cal_gyros.csv"
CALIBRATION_FILES_ACCELEROMETER = ["cal_acc_x.csv", "cal_acc_y.csv", "cal_acc_z.csv"]
DATA_FILE = 'rot90.csv'
SEPARATOR = ','
TIMEFRAME = 10

QUAT_ROT_INIT = Quaternion([1, 0, 0, 0]) #Init the quaternion to the identity rotation (q = 1 + 0i + 0j + 0k = 1)
QUAT_MAG_AIM = Quaternion([0, 20.8e-3, 0.0696e-3, 43.3e-3])
W_ARRAY_INIT = [0, 0, 0]
TIME_INIT = 0

def get_cal_data(filename, offsets = None):
	data = []
	with open(filename) as f:
		line = f.readline()
		while line:
			data.append(line.split(";"))
			line = f.readline()
	print(data[0], data[1])
	return [[float(r[i].replace(",", ".")) for r in data[1:]] for i in range(len(data[0]))]

def detect_sat(d, std_dev, thres=20):
	count = 0
	c = [0]
	for i in range(1, len(d)):
		if abs(d[i - 1] - d[i]) < std_dev/10:
			count += 1
			if count > thres:
				return True
		else:
			count = 0
	return False

def inv_sqrt(x):
	return 1/np.sqrt(x)

def norm(array):
	return sqrt(sum([x**2 for x in array]))

def normalize(array):
	return [x/norm(array) for x in array]

def get_raw_data(f):
	line, _ = f.recvfrom(8192)
	line = line.decode("ascii")

	if not line:
		return False

	data = [float(x.replace(",", ".")) for x in line.split(SEPARATOR)]

	return [
		data[0], # Time
		data[6:9], # Rotation velocity
		data[2:5], # Linear acceleration
		data[10:13]
	]

def jacobian(q, d):
	return np.matrix([	[0, 0, 0, 0],
				[2 * d[2] * q[3] - 2 * d[3] * q[2], 2 * d[2] * q[2] - 2 * d[3] * q[3], -4 * d[1] * q[2] + 2 * d[2] * q[1] - 2 * d[3] * q[0], -4 * d[0] * q[3] + 2 * d[2] * q[0] - 2 * d[3] * q[1]],
				[-2 * d[1] * q[3] + 2 * d[3] * q[1], 2 * d[1] * q[2] - 4 * d[2] * q[1] + 2 * d[3] * q[0], 2 * d[1] * q[2] + 2 * d[3] * q[1], -2 * d[1] * q[1] - 4 * d[2] * q[3] + 2 * d[3] * q[2]],
				[2 * d[1] * q[2] - 2 * d[2] * q[1], 2 * d[1] * q[3] - 2 * d[2] * q[0] - 4 * d[3] * q[1], 2 * d[1] * q[0] + 2 * d[2] * q[3] - 4 * d[3] * q[2], 2 * d[1] * q[1] + 2 * d[2] * q[2]]])

def quat_from_static(a, b, a_aim, b_aim, old, dt):
	f_a = old.inverted()*a_aim*old - a_aim
	f_b = old.inverted()*b_aim*old - b_aim

	j_a = jacobian(a, a_aim)
	j_b = jacobian(b, b_aim)


	f = np.bmat([[list(f_a)], [list(f_b)]]).transpose()
	j = np.bmat([[j_a.transpose(), j_b.transpose()]])

	delta_f = np.dot(j, f)
	return old - dt * Quaternion(delta_f).normalized()

print("Getting calibration data")

print("Gyros ---")
cal_data_gyros = get_cal_data(CALIBRATION_FILE_GYROS)
offsets_gyros = []

for data_set in cal_data_gyros[1:]:
	offsets_gyros.append(sum(data_set) / len(data_set))

print("Offsets: ", offsets_gyros)

print("Accelerometer ---")
offsets_accelerometer  = []
cal_values = [[], [], []]

for i, f in enumerate(CALIBRATION_FILES_ACCELEROMETER):
	cal_data_accelerometer = get_cal_data(f)
	for j in range(len(cal_values)):
		if j != i:
			cal_values[j] += cal_data_accelerometer[j]

offsets_accelerometer = [sum(c) / len(c) for c in cal_values]

print("Offsets: ", offsets_accelerometer)

print("Calibration done\n")

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM) # Setting up UDP Socket
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
s.bind(("", 5555))

plt.ion()
fig, ax = plt.subplots(2, sharex=True)
lines_rot = ax[0].plot([],[],[],[],[],[])
lines_acc = ax[1].plot([],[],[],[],[],[])

ax[0].legend(lines_rot, ["phi", "theta", "psi"])
ax[1].legend(lines_acc, ["x", "y", "z"])
plt.show()

ax[0].set_ylim(-pi, pi)
ax[1].set_ylim(-20, 20)

angles = [[], [], []]
acc = [[], [], []]
t = []
c = 0
previous_data = [(TIME_INIT, QUAT_ROT_INIT, W_ARRAY_INIT)]
terminate = False

q_static = QUAT_ROT_INIT

while not terminate:
	raw_data = get_raw_data(s)

	if not raw_data: # Escape the loop if there's no more data
		terminate = True
		continue

	for i in range(len(offsets_gyros)):
		raw_data[1][i] -= offsets_gyros[i]

	for i in range(len(offsets_accelerometer)):
		raw_data[2][i] -= offsets_accelerometer[i]

	t.append(raw_data[0])


	p_t, p_quat, p_w = previous_data[-1] # Get the old data

	new_quat_rot = Quaternion(normalize(raw_data[1]), norm(raw_data[1]) * (raw_data[0] - t[0] - p_t)) # Calculate the new quaternion

	q = p_quat * new_quat_rot

	q_static = quat_from_static(Quaternion([0] + raw_data[2]).normalized(), Quaternion([0, 0, 0, 1]), Quaternion([0] + raw_data[3]).normalized(), QUAT_MAG_AIM.normalized(), q_static, raw_data[0] - t[0] - p_t).normalized()


	previous_data.append((raw_data[0] - t[0], q, raw_data[1])) # Add the new data to the old

	phi, theta, psi = q_static.to_euler()
	angles[0].append(phi) # Compute Euler angles
	angles[1].append(theta)
	angles[2].append(psi)


	a = q_static * Quaternion([0] + raw_data[2]) * q_static.inverted() # Get acceleration in the fixed reference frame
	a[3] -= 9.813
	for i in range(len(acc)):
		acc[i].append(a[i+1])




# Display dara
	c += 1
	if c%5 == 0:
		for i in range(len(lines_rot)):
			lines_rot[i].set_data(t, angles[i])
		for i in range(len(lines_acc)):
			lines_acc[i].set_data(t, acc[i])
		plt.draw()
		ax[0].set_xlim(max(t[0], raw_data[0]-TIMEFRAME), raw_data[0])

plt.ioff()
plt.show()

