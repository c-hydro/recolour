
from converter import H16_Ind2ContTS as idx2cont

filefolder_in = '/share/VALIDATION_SM_HSAF/cell/tmp_idx/h16/'
filefolder_out = '/share/VALIDATION_SM_HSAF/cell/tmp_contiguous/h16/'

converter = idx2cont(filepath_in=filefolder_in, filepath_out=filefolder_out)
converter.conversion(cells=None)

