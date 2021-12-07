__all__ = ['contour_plot','plot_vectors']

from pylab import *

params =  {'text.usetex': True }
rcParams.update(params)



def get_ticks(values):
    vals = []

    if len(values) < 5:
        vals = values
    else:
        vals = []
        incr = (values[len(values)-1] - values[0])/4;

        for i in range(0,5):
            v = values[0] + incr * i
            if v < 10:
                vals.append(v)
            else:
                vals.append(int(v))

    return vals


def contour_plot(xs,ys,zs,lower,upper,step,title_text='',xtext='',ytext='',filename=''):
    figure()

    #x_row = map(lambda x: x[0],xs)
    #y_row = ys[0]
    extent = [0,10000,0,10000]
              
    im = imshow(zs, interpolation='bilinear', origin='lower',
                cmap=cm.gray,extent=extent)
    
    cset = contour(zs, arange(lower,upper,step),
                   origin='lower',
                   linewidths=2,extent=extent)
    
    clabel(cset,
           inline=1,
           fmt='%1.1f',
           fontsize=10)
    
    if xtext != '':
        xlabel(xtext)
        labels = get_ticks(x_row)
        yticks([0,2500,5000,7500,9200],tuple(labels))

        ylabel(ytext)
        labels = get_ticks(y_row)

        xticks([0,2500,5000,7500,9200],tuple(labels))

    else:
        axis('off')

    # hot()
    jet()
    colorbar()
    title(title_text)

    if filename != '':
        savefig(filename)
    else:
        show()


def plot_vectors(*vec):
    count = len(vec) * 100

    figure()

    indx = 1

    for v in vec:
        subplot(count + 10 + indx)
        indx += 1

        n  = v.n
        dt = v.dx

        t = arange(0,n)
        t = t * dt

        plot(t,v)
        # grid()

    show()
