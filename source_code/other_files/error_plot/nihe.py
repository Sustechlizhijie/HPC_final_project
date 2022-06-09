import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl
import warnings

x=[1e-5, 3e-5, 5e-5, 7e-5, 1e-4]
x1=np.arange(50,400,1)
y1=[0.00029782, 0.00074282, 0.00106975, 0.00131925, 0.00159996]
y2=np.array([395.6, 134.7, 105.7, 87.08, 52.91])
# y2=[395.6, 134.7, 105.7, 87.08, 52.91]

z = np.polyfit(y2, np.log(y1), 1)
p1=np.exp(z[0]*x1+z[1])
# p1= np.poly1d(z)
yvals = p1

plt.figure(num=1)
l1, = plt.plot(x1, yvals, color='black', linewidth=1.5, linestyle='--')
p1 = plt.scatter(y2, y1, s=30, alpha=1, color='black')
# leg1 = plt.legend(handles=[l1], labels=['analytical t1'],
#                         frameon=False, loc='best') #自动放在最好的位置,框不画出来
# plt.legend([p1, p2, p3], ['numerical t1'], loc='lower right',frameon=False, scatterpoints=1)
# plt.gca().add_artist(leg1)


#-----------------------
ax = plt.gca()  #获取坐标轴属性

#------------------------------------------设置坐标轴、坐标刻度参数
boxWidth = 1.5 #刻度的粗细
Lmajor = 5  #刻度整值上线的长度
Lminor = 3  #刻度分值上线的长度
xlabPad  = 10
ylabPad  = 10
# xlabel = r"$\mathdefault{u/u_0}$"
# ylabel = r"$\mathdefault{y/L}$"
xlabel = "time costing"
ylabel = "maximum error"
# xlimit = [50,400]
# ylimit = [0,0.0017]

# Hide the top and right spines of the axis--------选择边框是否可见
ax.spines['right'].set_visible(True)
ax.spines['top'].set_visible(True)
# Set the axis box line width---------------------设置坐标轴粗细
ax.spines['bottom'].set_linewidth(boxWidth)
ax.spines['left'].set_linewidth(boxWidth)
ax.spines['top'].set_linewidth(boxWidth)
ax.spines['right'].set_linewidth(boxWidth)
# Tick Parameters: Edit the major and minor ticks of the x and y axes---------坐标刻度属性
ax.xaxis.set_tick_params(which='major', size=Lmajor, width=boxWidth, direction='in',pad=xlabPad, top=False) #top边不用画出
ax.xaxis.set_tick_params(which='minor', size=Lminor, width=boxWidth, direction='in',pad=xlabPad, top=False)
ax.yaxis.set_tick_params(which='major', size=Lmajor, width=boxWidth, direction='in',pad=ylabPad, right=False)
ax.yaxis.set_tick_params(which='minor', size=Lminor, width=boxWidth, direction='in',pad=ylabPad, right=False)


plt.xlabel(xlabel)                 #坐标轴标签
plt.ylabel(ylabel)
# ax.set_xlim(xlimit[0], xlimit[1])  #坐标轴输出范围
# ax.set_ylim(ylimit[0], ylimit[1])

plt.show()