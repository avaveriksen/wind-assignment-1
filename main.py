from windtools import Tools, Interpolations # Our toolbox of helper functions

if __name__ == '__main__':
    '''main loop'''

    # Don't run this always, comment out when we have out tables
    #Interpolations.interp_tc() # Interpolate tables across tc

    #pn, pt, a, aprime, F = tools.bem_single_element(PROVIDE PARAMETERS)

    r=88.45
    c=1.18
    beta=-3.36
    V0= 4.0 # m/s
    omega=0.226116 # rad/s
    theta_p=-4.0 # degrees
    file = 'interpolated-tables\FFA_W3-0.241.csv'

    pn, pt, a, aprime, F = Tools.bem_single_element(r,c,beta,V0,omega,theta_p,file) # Run BEM on a single element, provide parameters within the function