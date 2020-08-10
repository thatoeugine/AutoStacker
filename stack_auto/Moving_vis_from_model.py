"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# This module places the WSCLean predicted visibilities from the Model column to Data column


"""
#--------------------------------------------------------------------------------------------



import casa
import numpy as np
import json
import numpy.ma as ma
import ConfigParser
import sys

if __name__=='__main__':
    config = ConfigParser.ConfigParser()
    config.read(sys.argv[-1])


    path = config.get('pipeline', 'data_path') 

    #Placing the predicted visibilities from the Model column to Data column
    casa.uvsub(path+config.get('simulator_params', 'msfile_name'), reverse = True)
    """
    tb.open(path+d["Simulator_params"]["msfile_name"].encode("ascii","ignore"),nomodify=False)  
    corrected_data=tb.getcol("CORRECTED_DATA")
    data = tb.getcol("DATA")
    adding_data = corrected_data + data
    tb.putcol("DATA", adding_data)
    #tb.putcol("CORRECTED_DATA", adding_data)

    tb.close()
    """


    def uvadd(vis):
        t = tbtool()
        casalog.post('Placing the predicted visibilities from the Model column to Data column', 'INFO')

        t.open(vis, nomodify=False)
        ram_restrict = 100000
        ranger = list(range(0,t.nrows(),ram_restrict))

        for colname in ['DATA']:
            if (colname in t.colnames()) and (t.iscelldefined(colname,0)):
                for j in ranger:
                    if j == ranger[-1]:
                        ram_restrict = t.nrows()%ram_restrict
                    a = t.getcol('DATA',startrow=j, nrow=ram_restrict, rowincr=1)
                    corrected_data = t.getcol('CORRECTED_DATA',startrow=j, nrow=ram_restrict, rowincr=1)
                    #model_mask = model == 0
                    #model = ma.array(data=model,mask=model_mask)
                    ### note that the use of np.isin cannot be used in CASA as it still uses numpy version v.1.11 whereas this came in v1.14
                    ### Had to use old horrible version (see line 51)
                    a = a+corrected_data
                    t.putcol('DATA',a,startrow=j, nrow=ram_restrict, rowincr=1)


    uvadd(path+config.get('simulator_params', 'msfile_name'))
    
