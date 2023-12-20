# python3 view.py
# sudo apt-get install pip
# sudo apt-get install python3-pip
# sudo pip3 install jupyter
# pip install h5py

import h5py
import os

file = open("output.txt", "w")

f = h5py.File("./MOF_PLT00000000.h5", "r")
count = 0

class Box:
    def __init__(self, F01, XCEN01, YCEN01, Order01, XSLOPE01, YSLOPE01, INTERCEPT01, F02, XCEN02, YCEN02, Order02, XSLOPE02, YSLOPE02, INTERCEPT02):
        self.F01 = F01
        self.XCEN01 = XCEN01
        self.YCEN01 = YCEN01
        self.YCEN02 = YCEN02
        self.Order01 = Order01
        self.XSLOPE01 = XSLOPE01
        self.YSLOPE01 = YSLOPE01
        self.INTERCEPT01 = INTERCEPT01
        self.F02 = F02
        self.XCEN02 = XCEN02
        self.Order02 = Order02
        self.XSLOPE02 = XSLOPE02
        self.YSLOPE02 = YSLOPE02
        self.INTERCEPT02 = INTERCEPT02
    
    def cout(self):
        file.write(str(self.F01)+ "\n")
        file.write(str(self.XCEN01)+ "\n")
        file.write(str(self.YCEN01)+ "\n")
        file.write(str(self.Order01)+ "\n")
        file.write(str(self.XSLOPE01)+ "\n")
        file.write(str(self.YSLOPE01)+ "\n")
        file.write(str(self.INTERCEPT01)+ "\n")
        file.write(str(self.F02)+ "\n")
        file.write(str(self.XCEN02)+ "\n")
        file.write(str(self.YCEN02)+ "\n")
        file.write(str(self.Order02)+ "\n")
        file.write(str(self.XSLOPE02)+ "\n")
        file.write(str(self.YSLOPE02)+ "\n")
        file.write(str(self.INTERCEPT02)+ "\n")

list = [
    "F01",
    "XCEN01",
    "YCEN01",
    "Order01",
    "XSLOPE01",
    "YSLOPE01",
    "INTERCEPT01",
    "F02",
    "XCEN02",
    "YCEN02",
    "Order02",
    "XSLOPE02",
    "YSLOPE02",
    "INTERCEPT02"
]

map = {}
cur = ""

for i_sets in f["Chombo_global"].attrs.values():
    print(i_sets)
    Chombo_spacedim=i_sets[0]

print("Chombo_spacedim=",Chombo_spacedim)

finest_level=0
for i_key in f.keys():
    finest_level=finest_level+1
    print(i_key)

finest_level=finest_level-2
print("finest_level=",finest_level)

for i in f["level_0"]["data:datatype=0"]:
    if (count % 576 == 0):
        cur = list[count // 576]
        map[cur] = []
        
    count += 1
    map[cur].append(i)

data = []

for i in range(576):
    # file.write("BOX #" + str(i) + "\n")
    data.append(Box(map[list[0]][i], map[list[1]][i], map[list[2]][i], map[list[3]][i], map[list[4]][i], map[list[5]][i], map[list[6]][i], map[list[7]][i], map[list[8]][i], map[list[9]][i], map[list[10]][i], map[list[11]][i], map[list[12]][i], map[list[13]][i]))
    # for att in list:
    #     file.write(str(att) + ": " + str(map[att][i]) + "\n")

for box in data:
    file.write("---------------------------------------------------------\n")
    box.cout()

file.close()
