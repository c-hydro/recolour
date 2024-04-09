import os
import pickle

folder_name = '/home/fabio/Desktop/PyCharm_Workspace/project/dte-ws/data_dynamic/ancillary/point_analysis_group/'
file_name_1 = 'analysis.soil_moisture.statistics_000305_000003.workspace'
file_name_2 = 'analysis.soil_moisture.statistics_000305_000006.workspace'

file_path_1 = os.path.join(folder_name, file_name_1)
file_obj_1 = pickle.load(open(file_path_1, "rb"))

file_path_2 = os.path.join(folder_name, file_name_2)
file_obj_2 = pickle.load(open(file_path_2, "rb"))

file_model_id_1 = file_obj_1['model_id']
file_model_id_2 = file_obj_2['model_id']


file_model_id_unique_1 = list(set(file_model_id_1))
file_model_id_unique_2 = list(set(file_model_id_2))

file_model_id_sum = file_model_id_unique_1 + file_model_id_unique_2

file_list = set(file_model_id_sum)

print('ciao')