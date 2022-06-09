import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# x = [1, 2, 4, 8]
# y = [8250, 9204, 21340, 16790]
x = [1, 2, 3, 4, 5]
y = [585.8, 3236, 9804, 17950, 12430]
yy = [585.8, 585.8*2, 585.8*3, 585.8*4, 585.8*5]


#-----------------------------------
plt.figure(num=1)
l1, = plt.plot(x, y, color='red', linewidth=1.5, linestyle='--')
p1 = plt.scatter(x, y, s=30, alpha=1, color='red')
l2, = plt.plot(x, yy, color='blue', linewidth=1.5, linestyle='solid')
p2 = plt.scatter(x, yy, s=30, alpha=1, color='blue')
leg1 = plt.legend(handles=[l1,l2], labels=['ture value', 'linear increase value'],
                        frameon=False, loc='best') #自动放在最好的位置,框不画出来
# plt.legend([p1, p2, p3], ['numerical t1'], loc='lower right',frameon=False, scatterpoints=1)
plt.gca().add_artist(leg1)
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
xlabel = "cores"
ylabel = "time(s)"
# xlimit = [-1e-3,1e-3]
# ylimit = [-1,1]

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