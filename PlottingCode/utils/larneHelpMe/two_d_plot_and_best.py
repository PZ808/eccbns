from overlap_calculator import overlap_calculator
from utils import contour_plot,find_best
from pan_et_al_values import *

from pylab import *

def two_d_plot_and_best(mass,approximant,order,params,target_dir):
    param = params[mass]
    eta   = param['eta']
    calc  = overlap_calculator('hybrid_merger',mass,approximant,order)

    # Calculate the grid
    # xs,ys,zs = calc.mass_freq_grid(m_steps=30,f_steps=30,eta=eta)

    #     xs,ys,zs = calc.mass_freq_grid(m_low  = 0,
    #                                    m_high = mass*2,
    #                                    m_steps=15,
    #                                    f_low  = param['f_high'] * (0.8),
    #                                    f_high = param['f_high'] * (1.2),
    #                                    f_steps=15,eta=eta)
    
    
    xs,ys,zs = calc.mass_eta_grid(param['f_high'],
                                  m_low  = 0,
                                  m_high = mass*2,
                                  m_steps=15,
                                  e_steps=15)
    

    # Find the highest overlap (z) and the mass (x) and f_c (y)
    # at which that overlap was found
    # TODO: Can I mark this point on the contour graph?
    best_x, best_y, best_z = find_best(xs,ys,zs)

    # Do the contour plot
    # title_text = 'Overlap between CC and ' + approximant + " PN " + str(order/2.0) + ", M = " + str(mass) + ", $\\eta=$" + str(eta)

    title_text = 'Overlap between CC and ' + approximant + " PN " + str(order/2.0) + ", M = " + str(mass) + ", $\f_c=$" + str(param['f_high'])

    filename = target_dir + '/CC_vs_' + approximant.replace(' ','_') + '_' + str(order/2.0).replace('.','_') + '_M' + str(mass) + '_eta.png'

    contour_plot(xs, ys, zs, 0.0, 1.0, 0.1, 
                 xtext='$f_c$',
                 ytext='Mass',
                 title_text=title_text,
                 filename=filename)


    #     # Now extract the f_c slice at best mass for a one-d plot...
    #     y_line = []
    #     z_line = []
    
    #     for i in range(0,len(xs)):
    #         for j in range(0,len(xs)):
    #             if xs[i][j] == best_x:
    #                 y_line.append(ys[i][j])
    #                 z_line.append(zs[i][j])
    
    
    #     # Or extract mass slice at best f_c...
    #     y_line = []
    #     z_line = []
    
    #     for i in range(0,len(xs)):
    #         for j in range(0,len(xs)):
    #             if ys[i][j] == best_y:
    #                 y_line.append(xs[i][j])
    #                 z_line.append(zs[i][j])
    
    #     # ... and plot it!
    #     figure()
    #     title(title_text)
    #     # xlabel('$f_c$')
    #     xlabel('Total Mass')
    #     ylabel('Overlap')
    #     plot(y_line,z_line)
    #     filename = target_dir + '/CC_vs_' + approximant.replace(' ','_') + '_' + str(order/2.0).replace('.','_') + '_M_' + str(mass) + '_m.png'
    #     savefig(filename)
    

    # Find the overlap using the values from Pan et. al.
    max = calc.calc_overlap_max(param['opt_mass'],
                                param['eta'],
                                param['f_high'])
    
    
    # Print a report line suitable to drop into the Wiki
    print "|", mass, "|", best_x, "|", best_y, "|", best_z, "|",param['opt_mass'],"|",param['eta'],"|",param['f_high'],"|",max,"|",param['overlap'],"|"


def print_header():
    print '^ Waveform M ^ M (found) ^ f<sub>c</sub> (found) ^ Overlap (found) ^ M (reported) ^ <html>&eta;</html> (reported) ^ f<sub>c</sub> (reported) ^ Overlap (calculated)  ^ Overlap (reported) ^'

def recreate_table_I():
    print_header()
    for mass in [10,20,30]:
	two_d_plot_and_best(mass,'Taylor F2',7,params_I,'table_I')
    

def recreate_table_II_35():
    print_header()
    for mass in [30,40,60,100]:
	two_d_plot_and_best(mass,'Taylor F2',7,params_II_35,'table_II_35')
    

def recreate_table_II_40():
    print_header()
    for mass in [30,40,60,100]:
	two_d_plot_and_best(mass,'Taylor F2',8,params_II_40,'table_II_40')



# recreate_table_I()
recreate_table_II_35()
# recreate_table_II_40()
