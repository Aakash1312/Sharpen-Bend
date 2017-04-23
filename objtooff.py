import sys
varr = []
farr = []
def getVertex(arr):
	t = (arr[1], arr[2], arr[3])
	return t
def getFace(arr):
	t = []
	size = len(arr)
	counter = 1
	while counter < size:
		splitarr = arr[counter].split('/')
		t.append(splitarr[0])
		counter = counter + 1
	return t
with open(sys.argv[1]) as f:
	while True:
		data = f.readline()
		if len(data) == 0:
			break
		arr = data.split()
		# print(arr)
		if len(arr) != 0 and arr[0] == "v":
			varr.append(getVertex(arr))
		if len(arr) != 0 and arr[0] == "f":
			farr.append(getFace(arr))

print("OFF")
print(len(varr), len(farr), "0")
for v in varr:
	print (float(v[0]),float(v[1]),float(v[2]))
for f in farr:
	print (len(f)," ",end='')
	for k in range(len(f)):
		if k == len(f) - 1:
			print(int(f[k]) - 1)
		else:
			print (int(f[k]) - 1," ",end='')
