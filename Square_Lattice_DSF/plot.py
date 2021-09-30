from parameters import *


def plot_fun(x_data,y_data,x_label,y_label,title,c):
    #if c=1, the plot is saved
    plt.plot(x_data,y_data)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    #plt.legend(["q = 0.2*pi","q = 1"])
    if(c==1):
        filename = "%s.csv" % title
        plt.savefig("%s.png",)
    plt.show()