import os
import re

import numpy
import matplotlib.pyplot as plt
D = 70
H = 10
v = 5
tf = 0.1
Dinter = 80
Dintra = 90
d = 115
stations = [[45.73, 45.26, 0.0], [1200, 700, 0], [-940, 1100, 0]]
sky_sites = [[-614, 1059, 24], [-934, 715, 12], [1073, 291, 37], [715, 129, 35], [186, 432, 21],  [-923, 632, 37], [833, 187, 24], [-63, 363, 11]]
def is_float(str1):
    if(str1.count('-')!=0 or str1.count('.')==0):
        return False
    return True
def distance(site1, site2):
    return ((site1[0] - site2[0])**2 + (site1[1] - site2[1])**2 + (site1[2] - site2[2])**2)**0.5

def is_time_not_legal(cost_time, dis):
    # print(cost_time, tf+ dis/10000)
    if(cost_time  < tf+ dis/10000.0):
        # print(cost_time, tf + dis / 10000)
        return True
    return False

def is_legal(path):
    last = path[0]
    for i in range(1, len(path)):
        now = path[i]
        # print(distance(now, last))
        if(distance(now, last) > d):
            return False
        last = now
    last = path[0]
    for i in range(1, len(path)):
        now = path[i]
        time_cost = now[3] - last[3]
        dis = distance(now, last)
        if(is_time_not_legal(time_cost, dis)):
            # print(len(path), i)
            # print(path[i])
            return False
        last = now
    return True

def result2(result_path):
    total_time = 0
    x = numpy.array([45.73, 1200, -940])
    y = numpy.array([45.26, 700, 1100])
    colors = ['r', 'b', 'y', 'g']
    plt.plot(x, y, 'r*', label='Ground Base Station')


    with open(result_path, 'r') as f:
        datas = f.readlines()
        station_end = []
        station_start = []
        count_wrong = 0
        time_cost_totall = 0
        pos = 0
        for index, line in enumerate(datas):

            if(index % 2 ==0):
                pos1 = int(line.split(',')[2])
                station_end = stations[pos1]
                pos2 = int(line.split(',')[1])
                station_start = stations[pos2]
                time = float(line.split(',')[3])
                time_cost_totall += time

                pass
            else:
                data = re.findall(r'[(](.*?)[)]', line)
                path = []
                # print('****')
                for site in data:
                    if (pos >= len(colors)):
                        pos = 0
                    site = numpy.array(site.split(',')).astype(numpy.float)
                    if(len(site)==3):
                        x, y = site[0]*v + site[1]*Dintra, site[2]*Dinter
                        plt.plot(x, y, colors[pos]+'o')
                        time = site[0]
                        path.append([x, y, H, site[0]])
                        # print(x, y, H)

                    if(len(site)==2):
                        x, y = sky_sites[int(site[1])][0], sky_sites[int(site[1])][1]
                        plt.plot(x, y, colors[pos] + '>')

                        time = site[0]
                        path.append([x, y, sky_sites[int(site[1])][2], float(site[0])])
                        # print(x, y, sky_sites[int(site[1])][2])
                pos += 1

                if(not is_legal(path)):
                    count_wrong+=1
                    # print("wrong")
                if(distance(station_end, path[-1]) > D):
                    print("wrong")
                if(distance(station_start, path[0]) > D):
                    print("wrong")
        print("wrong path: ", count_wrong)
        print("time cost:", time_cost_totall)
    x2 = numpy.array(range(0, 1200, 90))
    y2 = x2*0.567233+19.3204
    x3 = numpy.array(range(-940, 0, 90))
    y3 = x3*(-1.07001)+94.1915
    # print("时延: "+str(total_time))
    plt.plot(x2, y2, 'b--')
    plt.plot(x3, y3, 'b--')
    # # plt.savefig('result_062.png')
    plt.show()

path = './result.txt'
result2(path)