import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

data1 = np.loadtxt("explicit_CFL_0.1.dat", dtype=np.float)
data2 = np.loadtxt("implicit_CFL_1.0.dat", dtype=np.float)
y1=(data1[:])
y2=(data2[:])
#------------------------------算解析解
n=100
dx = 1.0/n
x = np.arange(0, 1+dx, dx)
y = np.sin(np.pi*x)/np.pi/np.pi

print(max(abs(y-y2)))
#
# #-----------------------------------
# plt.figure(num=1)
# l1, = plt.plot(x, y, color='black', linewidth=1.5, linestyle='solid')
# l2, = plt.plot(x, y1, color='red', linewidth=1.0, linestyle='--')
# l3, = plt.plot(x, y2, color='blue', linewidth=1.0, linestyle='-.')
# leg1 = plt.legend(handles=[l1,l2,l3], labels=['analytical solution','explicit method','implicit method'],
#                         frameon=False, loc='best') #自动放在最好的位置,框不画出来
#
# #-----------------------
# ax = plt.gca()  #获取坐标轴属性
#
# #------------------------------------------设置坐标轴、坐标刻度参数
# boxWidth = 1.5 #刻度的粗细
# Lmajor = 5  #刻度整值上线的长度
# Lminor = 3  #刻度分值上线的长度
# xlabPad  = 10
# ylabPad  = 10
# # xlabel = r"$\mathdefault{u/u_0}$"
# # ylabel = r"$\mathdefault{y/L}$"
# xlabel = r"$\mathdefault{X}$"
# ylabel = r"$\mathdefault{T}$"
# # xlimit = [-1e-3,1e-3]
# # ylimit = [-1,1]
#
# # Hide the top and right spines of the axis--------选择边框是否可见
# ax.spines['right'].set_visible(True)
# ax.spines['top'].set_visible(True)
# # Set the axis box line width---------------------设置坐标轴粗细
# ax.spines['bottom'].set_linewidth(boxWidth)
# ax.spines['left'].set_linewidth(boxWidth)
# ax.spines['top'].set_linewidth(boxWidth)
# ax.spines['right'].set_linewidth(boxWidth)
# # Tick Parameters: Edit the major and minor ticks of the x and y axes---------坐标刻度属性
# ax.xaxis.set_tick_params(which='major', size=Lmajor, width=boxWidth, direction='in',pad=xlabPad, top=False) #top边不用画出
# ax.xaxis.set_tick_params(which='minor', size=Lminor, width=boxWidth, direction='in',pad=xlabPad, top=False)
# ax.yaxis.set_tick_params(which='major', size=Lmajor, width=boxWidth, direction='in',pad=ylabPad, right=False)
# ax.yaxis.set_tick_params(which='minor', size=Lminor, width=boxWidth, direction='in',pad=ylabPad, right=False)
#
#
# plt.xlabel(xlabel)                 #坐标轴标签
# plt.ylabel(ylabel)
# # ax.set_xlim(xlimit[0], xlimit[1])  #坐标轴输出范围
# # ax.set_ylim(ylimit[0], ylimit[1])
#
# plt.show()