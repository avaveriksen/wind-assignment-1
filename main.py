from windtools import Tools, Interpolations # Our toolbox of helper functions

if __name__ == '__main__':
    '''main loop'''

    # Don't run this always, comment out when we have out tables
    Interpolations.interp_tc() # Interpolate tables across tc

    #pn, pt, a, aprime, F = tools.bem_single_element(PROVIDE PARAMETERS)

