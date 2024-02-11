from matplotlib import pyplot as plt
import numpy as n
f=open("leaky_1_25_si3n4_metalens_phase.txt","r")
scan=False
data= []
while True:
    line=f.readline()
    if not line:
        break
    if "(101,101)" in line:
        scan=True
        continue
    if line=="":
        continue
    if scan and len(line.split())>5:
        data.append([float(x) for x in line.split()])
    else:
        continue
data=n.array(data)
# import pdb; pdb.set_trace()
# # for i in range(int(data.shape[1]*0.1)):
# vec=data[:,24]
# n.savetxt("vec_polot.txt", vec)
# plt.plot(vec)

# # plt.imshow(data.T, cmap="hot")
# plt.legend()
# plt.show()
for i in range(data.shape[1]):
    vec=data[:,i]
    plt.plot(vec,label=str(i))

plt.legend()
plt.show()