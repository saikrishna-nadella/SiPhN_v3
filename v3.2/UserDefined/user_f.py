# -*- coding: utf-8 -*-

def f_udc(Re):
    """
      User defined correlation for friction factor in terms of eq. (1)

      f = p/(Re)^b    ---- (1)
      Data given in terms of arguments: 
      Re = Reynolds number of flow in section

      function returns p and b
    
    """
    # Add your correlation here -------------------------------------------
    p = 0.316
    b = 0.25

    # --------- User correlation ends -------------------------------------
    return (p,b)

