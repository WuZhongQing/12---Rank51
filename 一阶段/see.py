import matplotlib.pyplot as plt
import numpy
from matplotlib.font_manager import FontProperties

font = FontProperties(fname='/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf', size=14)
v = 5
Dinter = 80
Dintra = 90
def figure1():
    x = numpy.array([45.73, 1200, -940])
    y = numpy.array([45.26, 700, 1100])
    x1 = numpy.array(range(-1000, 1200, 90))
    print(len(x1))
    y1 = numpy.array(range(-800, 1200, 80))
    print(len(y1))
    plt.plot(x, y, 'r*')

    plt.plot(x1, y1, 'go')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(c='b')
    plt.savefig("grid.png")
    plt.show()
def design1():
    x = numpy.array([45.73, 1200])
    y = numpy.array([45.26, 700])
    x1 = numpy.array(range(0, 1200, 90))
    y1 = numpy.array(range(0, 1200, 80))
    plt.figure(figsize=(12, 7), dpi=90)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.plot(x, y, 'r*', label='Ground Base Station', linewidth=2, alpha=1)
    for i in y1[:-1]:
        plt.plot(x1, numpy.ones([len(x1)])*i, 'g>')
    plt.plot(x1, numpy.ones([len(x1)]) * y1[-1], 'g>', label='Unmanned Aerial Vehicle ',alpha=1)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.savefig('design1.png')
    plt.show()
def design2():
    x = numpy.array([-940, 45.73, 1200])
    y = numpy.array([1100, 45.26, 700])
    x1 = numpy.array(range(90, 1200, 90))
    y1 = numpy.array(range(0, 1200, 80))
    x2 = numpy.arange(-1000, 0, 90)
    y2 = numpy.arange(0, 1200, 80)
    plt.figure(figsize=(12, 7), dpi=90)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.plot(x, y, 'r*--', label='Ground Base Station', linewidth=2, alpha=1)
    for i in y1[:-1]:
        plt.plot(x1, numpy.ones([len(x1)])*i, 'g>')
    plt.plot(x1, numpy.ones([len(x1)]) * y1[-1], 'g>', label='Unmanned Aerial Vehicle ',alpha=1)
    for i in y2[:-1]:
        plt.plot(x2, numpy.ones([len(x2)])*i, 'g>')
    plt.plot(x2, numpy.ones([len(x2)]) * y2[-1], 'g>', label='Unmanned Aerial Vehicle ', alpha=1)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    # plt.savefig('design1.png')
    plt.show()
def design3():
    x = numpy.array([45.73, 1200])
    y = numpy.array([45.26, 700])
    x1 = numpy.array(range(0, 1200, 90))
    y1 = numpy.array(range(0, 1200, 80))
    plt.figure(figsize=(12, 7), dpi=90)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.plot(x, y, 'r*--', label='Ground Base Station', linewidth=2, alpha=1)
    for i in y1[:-1]:
        plt.plot(x1, numpy.ones([len(x1)])*i, 'g>')
    plt.plot(x1, numpy.ones([len(x1)]) * y1[-1], 'g>', label='Unmanned Aerial Vehicle ',alpha=1)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.savefig('design1.png')
    plt.show()
def result1(result_path):
    total_time = 0
    x = numpy.array([45.73, 1200, -940])
    y = numpy.array([45.26, 700, 1100])
    colors = ['r', 'b', 'y', 'g', 'g']
    plt.plot(x, y, 'r*', label='Ground Base Station')
    with open(result_path, 'r') as f:
        datas = f.readlines()
        datas = [data.strip().split(',') for data in datas]

        datas = [[data.replace('(', '').replace(')', '') for data in data_1] for data_1 in datas]
        for data in datas:
             if(len(data)==4):
                 data = numpy.array(data).astype(numpy.float)
                 total_time+= data[-1]
             else:
                 data = numpy.array(data).astype(numpy.double)
                 x_1 = []
                 y_1 = []
                 for i in range(0, len(data), 3):
                     t, m,n = data[i], data[i+1], data[i+2]
                     x1 = v*t+m*Dintra
                     y1 = n*Dinter
                     x_1.append(x1)
                     y_1.append(y1)

                     plt.plot(x1, y1, 'g>')
                     # x = numpy.arange(x1 - 125, x1 + 125, 0.01)
                     # y = y1 + numpy.sqrt(125 ** 2 - (x - y1) ** 2)
                     # plt.plot(x, y)
                     # plt.plot(x, -y)
                 # plt.plot(numpy.array(x_1), numpy.array(y_1), 'g>--')
    x2 = numpy.array(range(0, 1200, 90))
    y2 = x2*0.567233+19.3204
    x3 = numpy.array(range(-940, 0, 90))
    y3 = x3*(-1.07001)+94.1915
    print("时延: "+str(total_time))
    plt.plot(x2, y2, 'b--')
    plt.plot(x3, y3, 'b--')
    # plt.savefig('result_062.png')
    plt.show()

if __name__ == '__main__':
    result_path = "./result.txt"

    result1(result_path)
    # design2()